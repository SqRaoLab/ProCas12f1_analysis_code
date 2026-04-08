package main

import (
	"fmt"
	"strings"
)

// QualityScore 表示质量分数的类型和值
type QualityScore struct {
	Value   int  // ASCII转换后的质量分数（Phred分数）
	Char    rune // 原始ASCII字符
	IsValid bool // 是否有效
}

// QualityStats 质量统计信息
type QualityStats struct {
	MeanScore     float64     // 平均质量分数
	MinScore      int         // 最低质量分数
	MaxScore      int         // 最高质量分数
	Q20Count      int         // 质量≥20的碱基数
	Q30Count      int         // 质量≥30的碱基数
	Q20Percent    float64     // 质量≥20的碱基百分比
	Q30Percent    float64     // 质量≥30的碱基百分比
	LowQualityPos []int       // 低质量碱基的位置（可自定义阈值）
	ScoreHist     map[int]int // 质量分数分布直方图
}

// Phred33ToScore 将Phred+33编码的ASCII字符转换为质量分数
func Phred33ToScore(qual rune) (int, error) {
	score := int(qual) - 33
	if score < 0 || score > 93 { // Phred+33范围：0-93
		return 0, fmt.Errorf("无效的质量分数ASCII码: %c (%d)", qual, int(qual))
	}
	return score, nil
}

// Phred64ToScore 将Phred+64编码的ASCII字符转换为质量分数
func Phred64ToScore(qual rune) (int, error) {
	score := int(qual) - 64
	if score < 0 || score > 62 { // Phred+64范围：0-62
		return 0, fmt.Errorf("无效的质量分数ASCII码: %c (%d)", qual, int(qual))
	}
	return score, nil
}

// DetectPhredEncoding 检测质量字符串的编码格式（Phred+33 或 Phred+64）
func DetectPhredEncoding(qual string) string {
	if len(qual) == 0 {
		return "unknown"
	}

	// 检查ASCII范围
	minChar := 255
	maxChar := 0

	for _, c := range qual {
		if int(c) < minChar {
			minChar = int(c)
		}
		if int(c) > maxChar {
			maxChar = int(c)
		}
	}

	// Phred+33: ASCII 33-126 (通常33-74)
	// Phred+64: ASCII 64-126 (通常64-104)

	if minChar >= 33 && maxChar <= 74 {
		return "phred33" // 旧版本Illumina 1.8+ 使用Phred+33
	} else if minChar >= 64 && maxChar <= 126 {
		return "phred64" // 旧版本Illumina 1.5-1.7 使用Phred+64
	} else if minChar >= 33 && maxChar <= 126 {
		// 如果范围较广，看中位数
		medianChar := (minChar + maxChar) / 2
		if medianChar < 64 {
			return "phred33"
		}
		return "phred64"
	}

	return "unknown"
}

// Read 表示 FASTQ reads
type Read struct {
	ID, Seq, Qual string
}

func (r *Read) Id() string {
	id_ := strings.Split(r.ID, " ")[0]
	if strings.Contains(id_, "/1") || strings.Contains(id_, "/2") {
		id_ = strings.Split(id_, "/")[0]
	}
	return id_
}

func (r *Read) String() string {
	return fmt.Sprintf("%s\n%s\n+\n%s\n", r.ID, r.Seq, r.Qual)
}

// GetQualityScores 获取所有碱基的质量分数（自动检测编码）
func (r *Read) GetQualityScores() ([]QualityScore, error) {
	if len(r.Seq) != len(r.Qual) {
		return nil, fmt.Errorf("序列长度(%d)和质量字符串长度(%d)不匹配", len(r.Seq), len(r.Qual))
	}

	// 自动检测编码格式（Phred+33 或 Phred+64）
	encoding := DetectPhredEncoding(r.Qual)

	scores := make([]QualityScore, len(r.Qual))
	for i, c := range r.Qual {
		var score int
		var err error

		switch encoding {
		case "phred33":
			score, err = Phred33ToScore(c)
		case "phred64":
			score, err = Phred64ToScore(c)
		default:
			return nil, fmt.Errorf("无法检测质量分数编码")
		}

		if err != nil {
			scores[i] = QualityScore{Char: c, IsValid: false}
		} else {
			scores[i] = QualityScore{Value: score, Char: c, IsValid: true}
		}
	}

	return scores, nil
}

// CalculateQualityStats 计算质量统计信息
func (r *Read) CalculateQualityStats(thresholds ...int) (*QualityStats, error) {
	scores, err := r.GetQualityScores()
	if err != nil {
		return nil, err
	}

	// 设置低质量阈值（默认20）
	lowQualThreshold := 20
	if len(thresholds) > 0 {
		lowQualThreshold = thresholds[0]
	}

	stats := &QualityStats{
		MinScore:      100,
		MaxScore:      -1,
		ScoreHist:     make(map[int]int),
		LowQualityPos: []int{},
	}

	totalValid := 0
	sumScores := 0

	for i, sq := range scores {
		if !sq.IsValid {
			continue
		}

		totalValid++
		sumScores += sq.Value

		// 更新最大最小值
		if sq.Value < stats.MinScore {
			stats.MinScore = sq.Value
		}
		if sq.Value > stats.MaxScore {
			stats.MaxScore = sq.Value
		}

		// 直方图
		stats.ScoreHist[sq.Value]++

		// 低质量碱基统计
		if sq.Value >= 20 {
			stats.Q20Count++
		}
		if sq.Value >= 30 {
			stats.Q30Count++
		}
		if sq.Value < lowQualThreshold {
			stats.LowQualityPos = append(stats.LowQualityPos, i)
		}
	}

	if totalValid > 0 {
		stats.MeanScore = float64(sumScores) / float64(totalValid)
		stats.Q20Percent = float64(stats.Q20Count) / float64(totalValid) * 100
		stats.Q30Percent = float64(stats.Q30Count) / float64(totalValid) * 100
	}

	return stats, nil
}

// FilterByQuality 根据质量过滤reads
func (r *Read) FilterByQuality(minMeanQual int, minQ30Percent float64) (bool, error) {
	stats, err := r.CalculateQualityStats()
	if err != nil {
		return false, err
	}

	// 检查平均质量
	if stats.MeanScore < float64(minMeanQual) {
		return false, nil
	}

	// 检查Q30比例
	if stats.Q30Percent < minQ30Percent {
		return false, nil
	}

	return true, nil
}
