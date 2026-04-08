package main

import (
	"fmt"
	"strings"
)

// subsetString 采用坐标从bytes中截取对应的值
func subsetString(seq string, start, end int) string {
	if start >= 0 && end <= len(seq) && start < end {
		return seq[start:end]
	}
	//sugar.Debugf("subsetString start and end pos error: %d-%d", start, end)
	return ""
}

// SeqComponent 记录序列中的各种元件
type SeqComponent struct {
	PAM      string
	Edit     string
	Scaffold string
	Before   string
	Target   string
	After    string
	Spacer   string
	UMI      string
}

// Header used for output format
func (s *SeqComponent) Header(umi bool) []string {
	headers := []string{"scaffold", "spacer", "before", "pam", "target", "after", "edit"}
	if umi {
		headers = append(headers, "umi")
	}

	return headers
}

// String used for output results
func (s *SeqComponent) String(umi bool) []string {
	vals := []string{s.Scaffold, s.Spacer, s.Before, s.PAM, s.Target, s.After, s.Edit}

	if umi {
		vals = append(vals, s.UMI)
	}

	return vals
}

// findAllIndexes 返回子字符串在主字符串中所有出现的起始索引
func findAllIndexes(text, substr string) []int {
	var indexes []int
	start := 0
	for {
		if start >= len(text) {
			break
		}

		index := strings.Index(text[start:], substr)
		if index == -1 {
			break
		}
		// 记录全局位置
		globalIndex := start + index
		indexes = append(indexes, globalIndex)
		// 从下一个位置继续查找（避免重复匹配同一个起始点）
		start = start + index + 1
	}
	return indexes
}

// 解析测序序列
func (designSize *DesignSize) decodeString(seq, umi string, library *PoolLibrary) []*SeqComponent {

	data := make([]*SeqComponent, 0)
	res := &SeqComponent{}

	// 预计算 designSize 相关偏移
	if umi != "" {
		umiAnchors := strings.Split(umi, "|")
		indexes := make([]int, 0)
		for _, ua := range umiAnchors {
			indexes = append(indexes, strings.Index(seq, ua))
		}

		start, end := 0, 0
		if indexes[0] > 0 {
			start = indexes[0] + len(umiAnchors[0])
		}
		if indexes[1] > 0 {
			end = indexes[1]
		}

		if start > 0 && end > 0 {
			res.UMI = subsetString(seq, start, end)
		} else if start > 0 {
			res.UMI = subsetString(seq, start, start+designSize.UMILen)
		} else if end > 0 {
			res.UMI = subsetString(seq, end-designSize.UMILen, end)
		}
	}

	// if scaffold do not represent in sequence then return
	//requiredLength := designSize.SpacerDistance + designSize.BeforeLen + designSize.PAMLen + designSize.TargetLen + designSize.SpacerLen + len(library.Scaffold)
	for _, idx := range findAllIndexes(seq, library.Scaffold) {
		// find spacer
		start := idx + len(library.Scaffold)
		end := start + designSize.SpacerLen

		segment := subsetString(seq, start, end)

		// if spacer exists
		if beforeSeqs, ok := library.Spacer[segment]; ok {
			res.Scaffold = library.Scaffold
			res.Spacer = segment

			// find before
			start = end + designSize.SpacerDistance
			end = start + designSize.BeforeLen
			res.Before = subsetString(seq, start, end)

			// get pam
			start = end
			end = start + designSize.PAMLen
			res.PAM = subsetString(seq, start, end)

			// Quality check
			_, ok := beforeSeqs[res.Before]
			if !ok {
				sugar.Debugf("%s not found in Before map", res.Before)
				continue
			}
			_, ok = library.Before[res.Before][res.PAM]
			if !ok {
				sugar.Debugf("%s not found in PAM map of before %s", res.PAM, res.Before)
				sugar.Debugf("%s\n%s\n%s", seq, res.Before, res.PAM)
				continue
			}

			res.Target = library.Target[res.Spacer]
			// get edit sequence

			// before get edit sequence, identify the behind sequence first

			for behind := range library.Behind[fmt.Sprintf("%s%s%s", res.Before, res.PAM, res.Spacer)] {
				idx := strings.Index(seq[end:], behind)

				if idx > 0 {
					res.Edit = subsetString(seq, end, end+idx)
				} else {
					res.Edit = subsetString(seq, end, end+designSize.TargetLen)
				}

				temp := &SeqComponent{
					PAM:      res.PAM,
					Edit:     res.Edit,
					Scaffold: res.Scaffold,
					Before:   res.Before,
					Target:   res.Target,
					After:    behind,
					Spacer:   res.Spacer,
					UMI:      res.UMI,
				}

				data = append(data, temp)
			}
		}
	}

	return data
}

// ReadPair 代表一对FASTQ读段
type ReadPair struct {
	ID string
	R1 string
	R2 string
}
