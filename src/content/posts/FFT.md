---
title: 快速傅里叶变换（FFT）原理与实现
published: 2026-05-04
description: 快速傅里叶变换（FFT）原理与实现
tags: [FFT]
category: 算法
slug: fft
image: ./images/fft.png
---

<iframe width="100%" height="468" src="//player.bilibili.com/player.html?isOutside=true&aid=203048598&bvid=BV1za411F76U&cid=262835785&p=1" scrolling="no" border="0" frameborder="no" framespacing="0" title="快速傅里叶变换(FFT)——有史以来最巧妙的算法？"allowfullscreen="true"></iframe>

---

## 什么是FFT？
快速傅里叶变换（Fast Fourier Transform，FFT）是一种高效计算离散傅里叶变换（DFT）的算法，由Cooley和Tukey于1965年提出。它将DFT的时间复杂度从O(N²)降低到O(N log N)，使傅里叶分析在实时信号处理中变得可行。

---

### 核心思想
FFT利用**分治策略**和**旋转因子的周期性**，将大点数的DFT分解为小点数的DFT组合计算。最常见的基2-FFT算法要求信号长度N为2的整数次幂。

---

## FFT应用场景
1. **信号处理**：音频/图像滤波、频谱分析
2. **通信系统**：OFDM调制解调、信道估计
3. **医学成像**：MRI、CT扫描重建
4. **雷达系统**：目标检测与测距
5. **数据压缩**：MP3/JPEG等压缩算法

---

## 算法原理

---

### 离散傅里叶变换（DFT）
对于N点离散信号x[n]，其DFT定义为：
$$X[k] = \sum_{n=0}^{N-1} x[n] \cdot e^{-j2\pi kn/N} \quad k=0,1,...,N-1$$

---
### FFT分治过程
1. **按奇偶分解**：
   $$X[k] = \sum_{m=0}^{N/2-1} x[2m] \cdot W_N^{2mk} + W_N^k \sum_{m=0}^{N/2-1} x[2m+1] \cdot W_N^{2mk}$$
2. **旋转因子性质**：
   $$W_N^{2mk} = W_{N/2}^{mk}$$
3. **递归计算**：
   $$X[k] = G[k] + W_N^k H[k]$$
   其中G[k]是偶数点DFT，H[k]是奇数点DFT
---

## 代码实现分析

---


### 数据结构定义
```c
typedef struct {
    float real;   // 实部
    float imag;   // 虚部
} compx;
```

---

### FFT核心函数
```c
/**
 * @brief 快速傅立叶变换-FFT
 * 
 * @param xin  待变换数组
 * @param Num  数组元素个数
 * @param xout 存储FFT结果的数组
 */
void FFT(const compx *xin, int Num, compx *xout)
{
    int f, m, LH, nm, i, k, j, L;
    double p, ps; // 双精度 p: 指数，ps: 指数
    int le, B, ip; // le: 级数间隔，B: 级数间隔的点数，ip: 级数间隔的点序号
    float pi; // 圆周率
    compx w, t; // w: 旋转因子，t: 临时变量

    // 初始化输出数组，复制输入数组的内容
    for (i = 0; i < Num; ++i) {
        xout[i].real = xin[i].real;
        xout[i].imag = xin[i].imag;
    }

    LH = Num / 2;
    f = Num;
    for (m = 1; (f = f / 2) != 1; m++);

    for (L = m; L >= 1; L--) {
        le = pow(2, L);
        B = le / 2;
        pi = 3.1415;
        for (j = 0; j <= B - 1; j++) {
            p = pow(2, m - L) * j;
            ps = 2 * pi / Num * p;
            w.real = cos(ps);
            w.imag = -sin(ps);
            for (i = j; i <= Num - 1; i += le) {
                ip = i + B;
                t = xout[i]; // 临时变量保存xout[i]
                xout[i].real = t.real + xout[ip].real; // 蝶形运算
                xout[i].imag = t.imag + xout[ip].imag;
                xout[ip].real = t.real - xout[ip].real;
                xout[ip].imag = t.imag - xout[ip].imag;
                xout[ip] = EE(xout[ip], w); // 乘以旋转因子
            }
        }
    }

    // 变址运算
    nm = Num - 2;
    j = Num / 2;
    for (i = 1; i <= nm; i++) {
        if (i < j) {
            t = xout[j];
            xout[j] = xout[i];
            xout[i] = t;
        }
        k = LH;
        while (j >= k) {
            j = j - k;
            k = k / 2;
        }
        j = j + k;
    }
}
```

---

### 关键运算函数
```c
/**
 * @brief 复数乘法
 * 
 * @param b1 : 复数b1
 * @param b2 : 复数b2
 * @return struct compx 
 */
compx EE(compx b1, compx b2)
{
    compx b3;
    b3.real = b1.real*b2.real - b1.imag*b2.imag;
    b3.imag = b1.real*b2.imag + b1.imag*b2.real;
    return (b3);
}
```

---

## 测试信号生成与验证

---

### 测试信号配置
```c
#define SIGNAL_LENGTH 1024 // 信号长度，必须为 2 的幂
#define SAMPLING_RATE 1024  // 采样率，单位 Hz
#define SIGNAL_FREQ 74.6      // 信号频率，单位 Hz

// 生成测试信号
void generate_test_signal(compx *signal) 
{
	int i; 
    for (i = 0; i < SIGNAL_LENGTH; ++i) {
        signal[i].real = sin(2 * M_PI * SIGNAL_FREQ * i / SAMPLING_RATE);
        signal[i].imag = 0.0f;
    }
}
```

---

### FFT结果分析
```c
void test_fft() {
    int i = 0;
    compx test_signal[SIGNAL_LENGTH];
    compx fft_result[SIGNAL_LENGTH];
    
    generate_test_signal(test_signal);
    FFT(test_signal, SIGNAL_LENGTH, fft_result);
    
    // 打印结果
    printf("FFT Result:\n");
    for (i = 0; i < SIGNAL_LENGTH / 2; ++i) {
        float magnitude = sqrt(fft_result[i].real * fft_result[i].real 
                             + fft_result[i].imag * fft_result[i].imag);
        float freq = (float)i * SAMPLING_RATE / SIGNAL_LENGTH;
        printf("Index: %d, Frequency: %f Hz, Magnitude: %f\n", i, freq, magnitude);
    }
    
    // 检查最大值
    float max_magnitude = 0.0f;
    int max_index = 0;
    
    for (i = 0; i < SIGNAL_LENGTH / 2; ++i) {
        float magnitude = sqrt(fft_result[i].real * fft_result[i].real 
                             + fft_result[i].imag * fft_result[i].imag);
        if (magnitude > max_magnitude) {
            max_magnitude = magnitude;
            max_index = i;
        }
    }
    
    // 使用max_index计算频率
    float max_freq = (float)max_index * SAMPLING_RATE / SIGNAL_LENGTH;
    printf("\nPeak at Index: %d, Frequency: %f Hz, Magnitude: %f\n", 
           max_index, max_freq, max_magnitude);
}
```

---

## 常见问题及解决方案

---

### 1. 频率分辨率不足
**问题**：无法区分相近频率的信号  
**解决**：增加采样点数
```c
#define SIGNAL_LENGTH 4096  // 提高频率分辨率
```

---

### 2. 频谱泄漏
**问题**：非整数周期截断导致能量分散  
**解决**：使用窗函数
```c
// 汉宁窗应用
for (int i = 0; i < SIGNAL_LENGTH; i++) {
    float window = 0.5 - 0.5 * cos(2 * M_PI * i / (SIGNAL_LENGTH-1));
    signal[i].real *= window;
}
```

---

### 3. 幅值校准
**问题**：FFT结果需要缩放才能反映真实幅值  
**解决**：对结果进行归一化
```c
magnitude = sqrt(real*real + imag*imag) / (SIGNAL_LENGTH / 2);
// 直流分量特殊处理
if (i == 0) magnitude /= 2; 
```

---

### 4. 复数信号处理
**问题**：实信号FFT的对称性  
**现象**：频谱关于Nyquist频率对称  
**利用**：只需计算前半部分频谱（0 ~ N/2）

---

## 实际运行结果
```
Peak at Index: 75, Frequency: 75.000000 Hz, Magnitude: 386.443970
```

---

## FFT优化技巧
1. **查表法**：预计算旋转因子
2. **并行化**：利用SIMD指令加速蝶形运算
3. **混合基算法**：支持非2的幂长度
4. **迭代实现**：避免递归开销

---

## 总结
FFT是数字信号处理的基石算法，理解其原理和实现细节对于开发高性能信号处理系统至关重要。通过合理选择参数（采样率、点数）和应用窗函数等技术，可以获得精确的频谱分析结果。
