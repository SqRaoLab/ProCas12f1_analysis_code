package main

import (
	"bufio"
	"compress/gzip"
	"io"
	"os"
	"strings"
	"sync"
	"time"

	"github.com/k0kubun/go-ansi"

	"github.com/schollz/progressbar/v3"
)

// --- 文件读取和并发处理 ---
// progressReader 是一个包装器，用于在读取时更新进度条
type progressReader struct {
	reader io.Reader // 原始 Reader (例如 *os.File)
	bar    *progressbar.ProgressBar
}

// Read 实现 io.Reader 接口
func (pr *progressReader) Read(p []byte) (n int, err error) {
	// 调用原始 Reader 的 Read 方法
	n, err = pr.reader.Read(p)
	// 更新进度条
	if n > 0 {
		_ = pr.bar.Add(n)
	}
	return n, err
}

func progressBar(size int64, showBytes bool) *progressbar.ProgressBar {
	if size > 0 {
		return progressbar.NewOptions64(
			size,
			progressbar.OptionSetDescription("Reading"),       // 设置描述
			progressbar.OptionSetWriter(ansi.NewAnsiStderr()), // 使用 ansi 防止奇怪的换行等显示错误
			progressbar.OptionUseANSICodes(true),              // avoid progressbar downsize error
			progressbar.OptionEnableColorCodes(true),
			progressbar.OptionSetElapsedTime(true),
			progressbar.OptionSetPredictTime(true),
			progressbar.OptionSetWriter(os.Stderr), // 避免 stdout 冲突
			progressbar.OptionShowBytes(showBytes), // 显示已读/总字节数
			progressbar.OptionShowCount(),          // 显示计数（可选）
			progressbar.OptionFullWidth(),          // 全宽显示（可选）
			progressbar.OptionThrottle(1*time.Second),
		)
	}
	return progressbar.NewOptions(-1,
		progressbar.OptionSetDescription("Reading"),
		progressbar.OptionSetWriter(ansi.NewAnsiStderr()), // 使用 ansi 防止奇怪的换行等显示错误
		progressbar.OptionUseANSICodes(true),              // avoid progressbar downsize error
		progressbar.OptionEnableColorCodes(true),
		progressbar.OptionSetItsString("items"),
		progressbar.OptionSpinnerType(14),      // 旋转光标
		progressbar.OptionSetWriter(os.Stderr), // 避免 stdout 冲突
		progressbar.OptionSetItsString("files"),
		progressbar.OptionThrottle(1*time.Second),
	)
}

type ReadPair struct {
	R1     *Read
	R2     *Read
	Sample string
	Switch bool
}

// readFastqInChannel 从gzip压缩的FASTQ文件读取记录，并发送到channel
func readFastq(filename string, recordChan chan<- Read, wg *sync.WaitGroup, bar *progressbar.ProgressBar) {

	file, err := os.Open(filename)
	if err != nil {
		sugar.Fatalf("Error opening file %s: %v", filename, err)
	}

	var gzReader *gzip.Reader
	if bar != nil {
		// --- 4. 创建进度条包装器 ---
		progressReaderInstance := &progressReader{
			reader: file, // 包装原始文件 Reader
			bar:    bar,  // 关联进度条
		}
		gzReader, err = gzip.NewReader(progressReaderInstance)
		if err != nil {
			sugar.Fatalf("Error creating gzip reader for %s: %v", filename, err)
		}
	} else {
		gzReader, err = gzip.NewReader(file)
		if err != nil {
			sugar.Fatalf("Error creating gzip reader for %s: %v", filename, err)
		}
	}

	scanner := bufio.NewScanner(gzReader)
	buffer := make([]string, 0, 4) // 缓冲4行

	for scanner.Scan() {
		line := scanner.Text()
		buffer = append(buffer, line)

		if len(buffer) == 4 { // 一个完整的FASTQ记录
			// 发送记录的副本到channel
			recordCopy := make([]string, 4)
			copy(recordCopy, buffer)
			recordChan <- Read{
				ID:   recordCopy[0],
				Seq:  strings.TrimSpace(recordCopy[1]),
				Qual: strings.TrimSpace(recordCopy[3]),
			}

			buffer = buffer[:0] // 重置缓冲区
		}
	}

	if err := scanner.Err(); err != nil {
		sugar.Fatalf("Error reading from %s: %v", filename, err)
	}

	if gzReader != nil {
		_ = gzReader.Close()
	}

	if bar != nil {
		_ = bar.Finish()
	}

	_ = file.Close()

	// 处理文件末尾可能存在的不完整记录（理论上不应存在）
	if len(buffer) > 0 {
		sugar.Infof("Warning: Incomplete FASTQ record found at end of %s", filename)
	}

	wg.Done()
	close(recordChan) // 读取完成后关闭 channel
}

func readFastqInChannel(filename []string, recordChan chan<- *ReadPair, workWg *sync.WaitGroup, bar *progressbar.ProgressBar) {
	var wg sync.WaitGroup

	// 启动读取 goroutine
	wg.Add(2)

	if len(filename) < 2 {
		sugar.Fatalf("pair end reads required")
	}

	r1Path, r2Path := filename[0], filename[1]

	// 解析
	r1Ch := make(chan Read, 50000) // channel
	r2Ch := make(chan Read, 50000)

	go readFastq(r1Path, r1Ch, &wg, bar)
	go readFastq(r2Path, r2Ch, &wg, bar)

	var reads1, reads2 = make(map[string]Read), make(map[string]Read)

	for {
		s1, ok := <-r1Ch
		if !ok {
			break
		}
		s2, ok := <-r2Ch
		if !ok {
			break
		}

		if s1.Id() == s2.Id() {
			recordChan <- &ReadPair{
				R1: &s1,
				R2: &s2,
			}
		} else {
			// 如果r2对应的r1已经被收录了，取出来用。没有就记录r2
			if s, ok := reads1[s2.Id()]; ok {
				recordChan <- &ReadPair{
					R1: &s,
					R2: &s2,
				}

				reads1[s1.Id()] = s1
			} else {
				sugar.Debugf("%s", s2.Id())
				reads2[s2.Id()] = s2
			}
			// 如果r1对应的r2已经被收录了，取出来用。没有就记录r1
			if s, ok := reads2[s1.Id()]; ok {
				recordChan <- &ReadPair{
					R1: &s1,
					R2: &s,
				}

				reads2[s2.Id()] = s2
			} else {
				sugar.Debugf("%s", s1.Id())
				reads1[s1.Id()] = s1
			}

			if len(reads1) > 100 || len(reads2) > 100 {
				sugar.Fatalf("too much mismatch reads name in %s and %s", r1Path, r2Path)
			}
		}
	}

	wg.Wait()
	workWg.Done()
}
