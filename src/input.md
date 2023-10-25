# kmc++  参数说明及流程





## input.yaml



```
# System
par_ltc: "BCC"      
nx: 20				# 原子xyz扩胞比例
ny: 20
nz: 20
par_compB: 0.01     #BCDV的比例，小数
par_compC: 0.00
par_compD: 0.0
par_compV: 0.002
Mutiple: false		#是否开启多空位，开启需要设置par_compV比例，否则默认位1个空位
MutipleSize: false	#是否开启多类型，默认为 ABV 三种类型，开启可以有 ABCDV 五种类型
# System paralle
paralle: false       #是否启用并行
numThreads: 4		#利用的cpu数
parallelx: 4		#区域划分， 最小为4，且为偶数。 xyz 是xyz轴的划分。
parallely: 4
parallelz: 4
# Simulation time parameters
par_time: 100.0     # 总的时间
time_conf: 10.0     # output conf 
par_step: 5001		# 总步数
step_log: 500000	# output cluster
Outposcar: 1000		# 每隔多少步输出 POSCAR

# Read File
read_file: false		#从文件读取
filepath: "site.dat"	#文件读取路径
filepackage: "data"		#运行后文件保存路径
# Kinetic parameters 
par_temp: 673.0                   #tempk
par_beta: 8.617332478e-5		  # kb
par_dis_rec: 2.598                # recombination distance
par_muvA: 6.46e+12                # vcc mu and Em


# Sets the starting number of atoms
par_radius_start: 5			#从多少个原子开始才算一个团簇

# Energy parameters
par_eSPA: 0.65        # Fe
par_eSPB: 0.56		  # Cu
par_eSPC: 0.70        # Ni
par_eSPD: 1.03		  # Mn

#A
par_eSPA1A: -0.778    # A-A
par_eSPA2A: -0.389

par_eSPA1B: -0.609    # A-B
par_eSPA2B: -0.344

par_eSPA1V: -0.161    # A-V
par_eSPA2V: -0.161

par_eSPA1C: -0.821    
par_eSPA2C: -0.399

par_eSPA1D: -0.648    
par_eSPA2D: -0.364

#B
par_eSPB1B: -0.581    # B-B
par_eSPB2B: -0.389

par_eSPB1V: -0.103    # B-V
par_eSPB2V: -0.206

par_eSPB1C: -0.692   
par_eSPB2C: -0.344

par_eSPB1D: -0.519    
par_eSPB2D: -0.249

#C
par_eSPC1C: -0.793  
par_eSPC2C: -0.389

par_eSPC1D: -0.831   
par_eSPC2D: -0.464

par_eSPC1V: -0.234  
par_eSPC2V: -0.351

#D
par_eSPD1D: -0.438   
par_eSPD2D: -0.389

par_eSPD1V: -0.151    
par_eSPD2V: -0.206

#V
par_eSPV1V: 0.315     # V-V
par_eSPV2V: -0.214
```





这个文档中的参数设置用于描述系统和能量参数，以及模拟时间和并行计算等方面。以下是每个部分的说明：

### System

- ```yaml
  - par_ltc: "BCC" 表示晶格类型为体心立方(BCC)结构
  - nx, ny, nz: 分别表示系统在三个维度上的大小为40
  - par_compB, par_compC, par_compD, par_compV: 分别表示组分B、C、D和V的比例
  - Mutiple: true 表示启用多空位缺陷模拟，false则表示只有一个空位
  - MutipleSize:  表示是否使用多类型模拟，默认只有ABV三种类型，如果True则需要配置C、D比例。
  ```

  

### System paralle

- ```yaml
  - paralle: false 表示禁用并行计算
  - parallelx, parallely, parallelz: 分别表示在x、y和z方向上的并行计算数量，**2 * 2 * 2**就表示的是分为8个部分
  ```

  

### Simulation time parameters

- ```yaml
  - par_time: 100.0 表示总模拟时间为100.0 s
  - time_conf: 10.0 表示输出构型的时间间隔为10.0 s
  - par_step: 5000001 表示模拟的总步数为5000001         
  - step_log: 200000 表示输出中间参数 的步数间隔为200000
  - Outposcar: 5000000 表示输出POSCAR的步数
  ```

  

### Read File

- ```yaml
  - read_file: false 表示从文件读取坐标体系
  - filepath: "site.dat" 表示文件路径为"site.dat"
  - filepackage: "zax" 表示保存文件的文件夹路径
  ```

  

### Kinetic parameters

- ```yaml
  - par_temp: 773.0 表示温度为773.0K
  - par_beta: 8.617332478e-5 表示β常数为8.617332478e-5
  - par_dis_rec: 2.598 表示复合距离为2.598
  - par_muvA: 6.46e+12 表示vcc mu和Em为6.46e+12
  ```

  

### Sets the starting number of atoms

- ```yaml
  par_radius_start: 5 表示团簇的起始原子数为5
  ```

  

### Energy parameters

以下部分需要vasp计算的出：

- ```yaml
  - par_eSPA, par_eSPB, par_eSPC, par_eSPD: 分别表示跃迁能量参数
  - par_eSPA1A, par_eSPA2A, par_eSPA1B, par_eSPA2B, par_eSPA1V, par_eSPA2V, par_eSPA1C, par_eSPA2C, par_eSPA1D, par_eSPA2D: A部分的能量参数
  - par_eSPB1B, par_eSPB2B, par_eSPB1V, par_eSPB2V, par_eSPB1C, par_eSPB2C, par_eSPB1D, par_eSPB2D: B部分的能量参数
  - par_eSPC1C, par_eSPC2C, par_eSPC1D, par_eSPC2D, par_eSPC1V, par_eSPC2V: C部分的能量参数
  - par_eSPD1D, par_eSPD2D, par_eSPD1V, par_eSPD2V: D部分的能量参数
  ```




## **串行程序流程图**



​												![整体2](E:\KMC\KMCplus\KMCplus\整体2.png)







​													![整体3](E:\KMC\KMCplus\KMCplus\整体3.png)



​				

​														![整体4](E:\KMC\KMCplus\KMCplus\整体4.png)







​																![整体5](E:\KMC\KMCplus\KMCplus\整体5.png)



## 并行流程



​																	![整体6](E:\KMC\KMCplus\KMCplus\整体6.png)















