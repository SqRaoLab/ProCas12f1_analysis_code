package main

import (
	"strings"
	"sync"
)

func processPairedEndBarcodes(r1, r2 *Read, rule SampleBarcodeRule) (bool, bool) {
	// 计算反向互补序列（提前计算，避免重复）
	r1Rev := reverseComplement(r1.Seq, true)
	r2Rev := reverseComplement(r2.Seq, true)
	r1QualRev := reverseComplement(r1.Qual, false)
	r2QualRev := reverseComplement(r2.Qual, false)

	processed, reverse := false, false

	// 1. 检查正向匹配：R1含Up，R2含Down
	forwardMatch := strings.HasPrefix(r1.Seq, rule.Up) && strings.HasPrefix(r2.Seq, rule.Down)
	if forwardMatch {
		processed, reverse = true, false
		// 完全匹配，不需要处理
	}

	// 2. 检查反向匹配：R1含Down，R2含Up（barcode互换）
	reverseMatch := strings.HasPrefix(r1.Seq, rule.Down) && strings.HasPrefix(r2.Seq, rule.Up)
	if reverseMatch {
		// barcode互换，可以选择交换R1和R2
		// 这里只是记录警告，不修改序列
		processed, reverse = true, true
	}

	// 3. 检查正向匹配 + 反向互补：R1的反向互补含Up，R2的反向互补含Down
	forwardMatchRev := strings.HasPrefix(r1Rev, rule.Up) && strings.HasPrefix(r2Rev, rule.Down)
	if forwardMatchRev {
		// 双端都需要反向互补
		r1.Seq, r2.Seq = r1Rev, r2Rev
		r1.Qual, r2.Qual = r1QualRev, r2QualRev
		processed, reverse = true, false
	}

	// 4. 检查反向匹配 + 反向互补：R1的反向互补含Down，R2的反向互补含Up
	reverseMatchRev := strings.HasPrefix(r1Rev, rule.Down) && strings.HasPrefix(r2Rev, rule.Up)
	if reverseMatchRev {
		// 双端反向互补且barcode互换
		r1.Seq, r2.Seq = r1Rev, r2Rev
		r1.Qual, r2.Qual = r1QualRev, r2QualRev
		processed, reverse = true, true
	}

	// 5. 无法匹配任何模式
	return processed, reverse
}

// trimBarcodes trim掉barcode序列
func trimBarcodes(r1, r2 *Read, rule SampleBarcodeRule) {
	r1.Seq = r1.Seq[len(rule.Up):]
	r1.Qual = r1.Qual[len(rule.Up):]
	r2.Seq = r2.Seq[len(rule.Down):]
	r2.Qual = r2.Qual[len(rule.Down):]
}

// processFilePair 处理一对 R1/R2 FASTQ 文件
func processFilePair(
	readChan chan *ReadPair,
	outputChan map[string]chan *ReadPair,
	rules []SampleBarcodeRule,
	workWg *sync.WaitGroup,
	minMeanQual int, minQ30Perc float32,
	trim bool,
) {
	for {
		readPair, ok := <-readChan
		if !ok {
			break
		}

		r1, r2 := readPair.R1, readPair.R2

		if ok, _ := r1.FilterByQuality(minMeanQual, float64(minQ30Perc)); !ok {
			continue
		}

		if ok, _ := r2.FilterByQuality(minMeanQual, float64(minQ30Perc)); !ok {
			continue
		}

		// iter barcodes for different cas protein
		for _, rule := range rules {

			processed, reverse := processPairedEndBarcodes(r1, r2, rule)
			if processed {
				if trim {
					trimBarcodes(readPair.R1, readPair.R2, rule)
				}
				readPair.Switch = reverse
				readPair.Sample = rule.Sample
				outputChan[rule.Sample] <- readPair
				break
			}
		}
	}

	workWg.Done()
}
