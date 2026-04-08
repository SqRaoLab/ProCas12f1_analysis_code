package main

import (
	"encoding/json"
	"os"
)

// SampleBarcodeRule 每个样本的 barcode 规则
type SampleBarcodeRule struct {
	Sample string
	Data   string
	Up     string
	Down   string
}

// loadBarcodes 从 JSON 文件加载 barcode 规则; key: sample; value of barcode rules
func loadBarcodes(path string) (map[string][]SampleBarcodeRule, map[string]string, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, nil, err
	}

	var raw map[string]map[string][][]string
	if err := json.Unmarshal(data, &raw); err != nil {
		return nil, nil, err
	}

	samples := make(map[string]string)
	rules := make(map[string][]SampleBarcodeRule)
	for sample, dataMap := range raw {
		for dataKey, pairs := range dataMap {
			for _, pair := range pairs {
				if len(pair) >= 2 {
					if _, ok := rules[dataKey]; !ok {
						rules[dataKey] = make([]SampleBarcodeRule, 0)
					}
					rules[dataKey] = append(rules[dataKey], SampleBarcodeRule{
						Sample: sample,
						Data:   dataKey,
						Up:     pair[0],
						Down:   pair[1],
					})
					samples[sample] = ""
				}
			}
		}
	}
	return rules, samples, nil
}
