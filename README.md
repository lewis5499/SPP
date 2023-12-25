# SPP

## 1.程序编译与运行

### 1.1 源码

项目采用CMake管理, CMakeLists.txt保存了项目配置信息;
./src目录下提供了所有源码

### 1.2 依赖库

除了基本的C++标准库之外, SPP项目无其他依赖库

### 1.3 运行结果

SPP定位结果保存至./log/
./log/中含有matlab代码, 用于轨迹可视化与rms/std可视化。

## 2 数据集

### 2.1 测试数据

./dataset/目录下提供了两份数据, "kinematic.gps"为动态的接收机数据，"202310301810.oem719"为静态的接收机数据

### 2.2 SPP自采集数据

如果以实时SPP模式解算，程序默认从网口采集播发的原始二进制数据，并命名为"raw.oem719"保存至./log/中
