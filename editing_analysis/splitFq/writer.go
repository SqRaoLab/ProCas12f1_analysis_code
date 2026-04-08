package main

import (
	"compress/gzip"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"sync"
)

// FileWriterCache 缓存输出文件句柄
type FileWriterCache struct {
	writers map[string]*gzip.Writer
	files   map[string]*os.File
}

func NewFileWriterCache() *FileWriterCache {
	return &FileWriterCache{
		writers: make(map[string]*gzip.Writer),
		files:   make(map[string]*os.File),
	}
}

func (c *FileWriterCache) GetWriter(sample, readID, outputDir string) (*gzip.Writer, error) {
	key := fmt.Sprintf("%s_%s.fq.gz", sample, readID)

	if readID == "" {
		key = fmt.Sprintf("%s.fq.gz", sample)
	}

	filename := filepath.Join(outputDir, key)
	if w, exists := c.writers[key]; exists {
		return w, nil
	}
	file, err := os.Create(filename)
	if err != nil {
		return nil, err
	}

	gz := gzip.NewWriter(file)
	c.writers[key] = gz
	c.files[key] = file
	return gz, nil
}

func (c *FileWriterCache) CloseAll() {
	for key, gz := range c.writers {
		_ = gz.Close()
		_ = c.files[key].Close()
	}
}

// reverseComplement 返回 DNA 序列的反向互补链
func reverseComplement(seq string, doComplement bool) string {
	// 定义互补映射
	complement := map[rune]rune{
		'A': 'T',
		'T': 'A',
		'G': 'C',
		'C': 'G',
		'a': 't', // 可选：支持小写
		't': 'a',
		'g': 'c',
		'c': 'g',
	}

	var builder strings.Builder
	builder.Grow(len(seq)) // 预分配空间，提升性能

	// 从后往前遍历，实现反向 + 互补
	for i := len(seq) - 1; i >= 0; i-- {
		r := rune(seq[i])

		if doComplement {
			if comp, exists := complement[r]; exists {
				builder.WriteRune(comp)
			} else {
				// 可选：处理非法字符，如 'N' 或其他
				builder.WriteRune(r) // 或 builder.WriteRune('N') 表示未知
			}
		} else {
			builder.WriteRune(r)
		}

	}

	return builder.String()
}

func writer(outChan chan *ReadPair, output string, workerWg *sync.WaitGroup, merged bool) {
	writerCache := NewFileWriterCache()

	count := 0
	for {
		out, ok := <-outChan
		if !ok {
			break
		}

		r1, r2 := out.R1, out.R2

		if out.Switch {
			r1, r2 = out.R2, out.R1
		}

		if merged {
			// merge PE into single fq file

			w1, err := writerCache.GetWriter(out.Sample, "", output)
			if err != nil {
				sugar.Fatal(err)
			}

			if strings.Contains(r1.ID, "E250040455L1C010R03400035384") {
				sugar.Infof(out.Sample)
			}

			_, _ = w1.Write([]byte(fmt.Sprintf(
				"%s\n%s%s\n+\n%s%s\n",
				r1.ID,
				r1.Seq, reverseComplement(r2.Seq, true),
				r1.Qual, reverseComplement(r2.Qual, false),
			)))

			if count%10000 == 0 {
				_ = w1.Flush()
			}
		} else {
			w1, err := writerCache.GetWriter(out.Sample, "R1", output)
			if err != nil {
				sugar.Fatal(err)
			}

			_, _ = w1.Write([]byte(fmt.Sprintf("%s\n%s\n+\n%s\n", r1.ID, r1.Seq, r1.Qual)))

			w2, err := writerCache.GetWriter(out.Sample, "R2", output)
			if err != nil {
				sugar.Fatal(err)
			}
			_, _ = w2.Write([]byte(fmt.Sprintf("%s\n%s\n+\n%s\n", r2.ID, r2.Seq, r2.Qual)))

			if count%10000 == 0 {
				_ = w1.Flush()
				_ = w2.Flush()
			}
		}

		count++
	}

	writerCache.CloseAll()
	workerWg.Done()
}
