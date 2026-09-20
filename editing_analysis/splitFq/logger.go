package main

import (
	"os"
	"time"

	"go.uber.org/zap"
	"go.uber.org/zap/zapcore"
	"gopkg.in/natefinch/lumberjack.v2"
)

var (
	sugar *zap.SugaredLogger
)

// NewEncoderConfig create EncoderConfig for zap
func newEncoderConfig() zapcore.EncoderConfig {
	return zapcore.EncoderConfig{
		// Keys can be anything except the empty string.
		TimeKey:        "T",
		LevelKey:       "L",
		NameKey:        "N",
		CallerKey:      "C",
		MessageKey:     "M",
		StacktraceKey:  "S",
		LineEnding:     zapcore.DefaultLineEnding,
		EncodeLevel:    zapcore.CapitalColorLevelEncoder,
		EncodeTime:     timeEncoder,
		EncodeDuration: zapcore.StringDurationEncoder,
		EncodeCaller:   zapcore.ShortCallerEncoder,
	}
}

// TimeEncoder format logger time format
func timeEncoder(t time.Time, enc zapcore.PrimitiveArrayEncoder) {
	enc.AppendString(t.Format("2006-01-02 15:04:05.000"))
}

// setLogger init logger
func setLogger(path string) {
	encoder := newEncoderConfig()

	consoleCore := zapcore.NewCore(
		zapcore.NewConsoleEncoder(encoder),
		zapcore.AddSync(os.Stderr),
		zap.InfoLevel,
	)

	// --- 2. 配置文件输出 (JSON Encoder) ---
	// 使用 lumberjack 实现日志文件轮转
	fileWriter := zapcore.AddSync(&lumberjack.Logger{
		Filename:   path, // 日志文件路径
		MaxSize:    10,   // 每个日志文件最大大小 (MB)
		MaxBackups: 5,    // 保留旧文件的最大个数
		MaxAge:     30,   // 保留旧文件的最大天数
		Compress:   true, // 是否压缩旧文件
	})

	// 文件日志通常使用 JSON 格式，便于解析
	fileEncoder := zapcore.NewConsoleEncoder(encoder)

	// 创建文件 Core，记录 DebugLevel 及以上级别的日志
	fileCore := zapcore.NewCore(fileEncoder, fileWriter, zapcore.DebugLevel)

	// --- 3. 组合 Cores ---
	// 使用 zapcore.NewTee 将多个 Core 组合起来
	// 这样日志会同时发送到所有匹配级别的 Core
	teeCore := zapcore.NewTee(consoleCore, fileCore)

	// 创建 Logger
	logger := zap.New(teeCore,
		// 添加调用者信息 (文件名和行号)
		zap.AddCaller(),
		// 添加堆栈跟踪 (对于 WarnLevel 及以上级别)
		zap.AddStacktrace(zapcore.WarnLevel),
	)

	defer logger.Sync()
	sugar = logger.Sugar()
}
