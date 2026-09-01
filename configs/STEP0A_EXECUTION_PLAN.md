# STEP0A Execution Plan — Draft, Not Authorized

本文件只定义 M00 的边界，不代表 pilot 已运行。

## 固定输入

四个 run 见 `02_protocol_and_design/input_manifests/STEP0A_INPUT_MANIFEST.csv`。不得用运行结果更好的其他 run 替换。

## 预运行仍需补齐

1. 每个下载对象的本地路径、文件大小和 SHA-256。
2. 采用的 read subset 或全量策略及随机性规则。
3. 人与鼠分类器/参考的名称、版本、下载来源和 SHA-256。
4. 工具版本、精确命令数组、阈值、线程数、资源预算和 watchdog 参数。
5. 建立真实的运行入口和代码哈希，以取代当前仅用于结构验证的执行计划哈希，并重新验证合同。

## 运行证据

每个 run 必须产生：下载/转换日志、输入哈希、工具/参考身份、命令、run-level 统计、分类、警告和退出码。四个 run 汇总后形成唯一 gate report。

## 禁止事项

- 不做 DEG、富集、WGCNA、免疫浸润、机器学习或图表故事化。
- 不把 metadata 标签当作 read-level PASS。
- 不把 inconclusive 当作 human PASS。
- 不从失败目录复制部分结果到 locked results。
