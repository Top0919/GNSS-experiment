# 引用包
import math
import numpy

# 函数
# 是否已经存在，用于判断后续冗余数据块
def isin(value1, value2, N_data):
    for i in range(len(N_data)):
        if N_data[i]['PRN号'] == value1:
            if N_data[i]['历元'] == value2:
                return True
    return False

# 返回头信息所占行数
def head_counts(lines_N):
    flag = 0
    for i in range(len(lines_N)):
        if lines_N[i].find('END OF HEADER') != -1:  # N文件中头信息终止标志
            flag = i + 1
    return flag

# 迭代进行偏近点角计算
def Ek(M, e):
    E = 0
    while True:
        e0 = E
        E = M + e * math.sin(e0)
        if abs(E - e0) < 0.00000001:  # 迭代结束条件
            break
    return E

# # cal2jd 将公历年月日时分秒转换到简化儒略日。
# def cal2mjd(calendar):
#     # 输入公历时间列表，返回儒略日
#     if len(calendar) < 6:
#         for i in range(len(calendar), 6):
#             calendar.append(0)
#     year = calendar[0]
#     month = calendar[1]
#     day = calendar[2] + (calendar[3] * 3600 + calendar[4] * 60 + calendar[5]) / 86400
#     y = year + 4800
#     m = month
#     if year < 0:
#         print('Year is wrong')
#         return False
#     if m <= 2:
#         # 1，2月视为前一年13，14月
#         m = m + 12
#         y = y - 1
#     e = math.floor(30.6 * (m + 1))
#     a = math.floor(y / 100)
#     # 教皇格雷戈里十三世于1582年2月24日以教皇训令颁布，将1582年10月5日至14抹掉。1582年10月4日过完后第二天是10月15日
#     if (year < 1582) or (year == 1582 and month < 10) or (year == 1582 and month == 10 and day < 15):
#         b = -38
#     else:
#         b = math.floor((a / 4) - a)
#     c = math.floor(365.25 * y)
#     jd = b + c + e + day - 32167.5
#     mjd = jd - 2400000.5
#     return mjd
#
# # cal2gps 将公历GPS时间转换到GPS周和周内的秒
# def cal2gps(calendar):
#     # cal2gps 将公历GPS时间转换到GPS周和周内的秒
#     mjd = cal2mjd(calendar)
#     # GPS时从1980年，MJD44244开始
#     e = mjd - 44244
#     week = math.floor(e / 7)    #周数
#     e = e - week * 7
#     return [week, round(e * 86400)]  # 返回列表，周和周内秒

# 0、读取N文件
def Read_N_file(path):
    with open(path, 'r') as f:  # 制度方式打开N文件
        lines_N = f.readlines()     # 逐行读取
    head_count = head_counts(lines_N)  # 计算头文件行数
    N_data = []  # 存储各个数据块数据
    Data_count = int((len(lines_N) - head_count) / 8)    # 数据块数量（去除头文件后）
    # 广播星历关键字,第一行单独定义方便循环
    N_keys = ['IODE', 'Crs', 'delta_N', 'M0',
              'Cuc', 'e', 'Cus', 'sqrt_a',
              't0', 'Cic', 'Ω0', 'Cis',
              'i0', 'Crc', 'w', 'Ωdot',
              'idot', 'l2上的码', 'gps_week', 'l2p码数据标记',
              'svacc', 'svhlth', 'TGD', 'IODC',
              'ttm', 'fi', 'spare1', 'spare2']
    for i in range(Data_count):  # 第一层不同数据块之间循环
        N_data_item = {}
        flag = 0  # 每行数据下标，0,1,2,3,
        for j in range(8):  # 每个数据块内的循环
            data_content = lines_N[head_count + 8 * i + j]  # 第i块j行数据内容
            N_data_item['数据号'] = i + 1
            if j == 0:  # line 1
                N_data_item['PRN号'] = data_content.strip('\n')[0:3].strip(' ')  # 左闭右开
                N_data_item['历元'] = data_content.strip('\n')[4:23].strip(' ')
                if isin(N_data_item['PRN号'], N_data_item['历元'], N_data):  # 消除冗余
                    break
                # 示例G03 2022 12  1 10  0  0 -3.74834053218D-04 -1.13686837722D-12  0.00000000000D+00
                # 科学计数法表示
                N_data_item['卫星钟偏差'] = str(float(data_content.strip('\n')[24:42][:-4]) * numpy.power(10.0, -int(data_content.strip('\n')[24:42][-2:])))
                N_data_item['卫星钟漂移'] = str(float(data_content.strip('\n')[43:61][:-4]) * numpy.power(10.0, -int(data_content.strip('\n')[43:61][-2:])))
                N_data_item['卫星钟漂移速度'] = str(float(data_content.strip('\n')[63:80][:-4]) * numpy.power(10.0,int(data_content.strip('\n')[63:80][-2:])))
            # 其他7行数据
            else:  # 第一数据
                if data_content.strip('\n')[6:23][-3] == '-':  # 科学计数法符号（负号）
                    N_data_item[N_keys[flag]] = str(float(data_content.strip('\n')[6:23][:-4]) * numpy.power(10.0, -int(data_content.strip('\n')[6:23][-2:])))
                else:
                    N_data_item[N_keys[flag]] = str(float(data_content.strip('\n')[6:23][:-4]) * numpy.power(10.0,int(data_content.strip('\n')[6:23][-2:])))
                    # 第二数据
                if data_content.strip('\n')[24:42][-3] == '-':
                    N_data_item[N_keys[flag + 1]] = str(float(data_content.strip('\n')[24:42][:-4]) * numpy.power(10.0,-int(data_content.strip('\n')[24:42][-2:])))
                else:
                    N_data_item[N_keys[flag + 1]] = str(float(data_content.strip('\n')[24:42][:-4]) * numpy.power(10.0,int(data_content.strip('\n')[24:42][-2:])))
                    # 第三数据
                if data_content.strip('\n')[43:61][-3] == '-':
                    N_data_item[N_keys[flag + 2]] = str(float(data_content.strip('\n')[43:61][:-4]) * numpy.power(10.0,-int(data_content.strip('\n')[43:61][-2:])))
                else:
                    N_data_item[N_keys[flag + 2]] = str(float(data_content.strip('\n')[43:61][:-4]) * numpy.power(10.0,int(data_content.strip('\n')[43:61][-2:])))
                    # 第四数据
                if data_content.strip('\n')[62:80][-3] == '-':
                    N_data_item[N_keys[flag + 3]] = str(float(data_content.strip('\n')[62:80][:-4]) * numpy.power(10.0,-int(data_content.strip('\n')[62:80][-2:])))
                else:
                    N_data_item[N_keys[flag + 3]] = str(float(data_content.strip('\n')[62:80][:-4]) * numpy.power(10.0,int(data_content.strip('\n')[62:80][-2:])))
                flag = flag + 4
        if len(N_data_item) > 4:
            N_data.append(N_data_item)
    return N_data


# 参数
# 2022年12月1日8时30分0秒 周四
t1 = ((4 * 24 + 8) * 60 + 30) * 60  # 周内秒
u = 3.986005E+14  # 地球引力常数
we = 7.2921150E-5  # 地球旋转速率rad/s
pi = 3.1415926535898  # 圆周率π
cal = [2022, 12, 1, 8, 30, 0]   # 2022年12月1日8点30日期

# 1、广播星历计算卫星坐标
def position(n_data):
    coord = []  # 用于存储每颗卫星坐标
    for i in range(len(n_data)):    # 依次计算每颗卫星坐标
        a = numpy.power(float(n_data[i]['sqrt_a']), 2)  # 长半轴半径a
        # 1计算归化时间,2022/12/1 8:30:00:00周四
        t = t1  # 用户时间
        toe = float(n_data[i]['t0'])    # 参考时间，从广播星历中获取
        tk = t - toe    # 归化时间
        # 2计算卫星平均角速度
        n0 = numpy.sqrt(u / numpy.power(a, 3))
        n = n0 + float(n_data[i]['delta_N'])
        # 3计算信号发射时刻平近点角
        m0 = float(n_data[i]['M0'])
        Ms = m0 + n * tk
        # 4计算信号发射时刻偏近点角(rad) Es=Ms+e*sin(Es)进行迭代
        e = float(n_data[i]['e'])
        Es = Ek(Ms, e)
        # 5真近点角
        fs = 2 * math.atan(numpy.sqrt((e + 1) / (1 - e)) * math.tan(Es / 2))    # 公式
        # 6升交点角距theta
        theta0 = float(n_data[i]['w']) + fs  # 广播星历文件数据
        # 7摄动改正项
        Cuc = float(n_data[i]['Cuc'])   # 从广播星历获取
        Cus = float(n_data[i]['Cus'])
        Crs = float(n_data[i]['Crs'])
        Cis = float(n_data[i]['Cis'])
        Crc = float(n_data[i]['Crc'])
        Cic = float(n_data[i]['Cic'])
        g_u = Cuc * math.cos(2 * theta0) + Cus * math.sin(2 * theta0)   # du=Cus*sin2θ+Cuc*2θ
        g_r = Crc * math.cos(2 * theta0) + Crs * math.sin(2 * theta0)   # dr=Crs*sin2θ+Crc*2θ
        g_i = Cic * math.cos(2 * theta0) + Cis * math.sin(2 * theta0)   # di=Cis*sin2θ+Cic*2θ
        # 8升交点角距、卫星的地心距离、轨道倾角
        theta = theta0 + g_u
        r = g_r + a * (1 - e * math.cos(Es))
        i1 = g_i + float(n_data[i]['i0']) + tk * float(n_data[i]['idot'])
        # 9计算信号发射时刻的升交点赤经
        Ω0 = float(n_data[i]['Ω0'])  # 广播星历文件数据
        Ωdot = float(n_data[i]['Ωdot'])
        Lk = Ω0 + (Ωdot - we) * tk - we * toe   # 升交点赤径，依据公式
        # 10计算卫星在瞬时地球坐标系中的坐标
        m1 = numpy.array([[math.cos(Lk), -math.sin(Lk) * math.cos(i1), math.sin(Lk) * math.sin(i1)],    # 矩阵1
                          [math.sin(Lk), math.cos(Lk) * math.cos(i1), -math.cos(Lk) * math.sin(i1)],
                          [0, math.sin(i1), math.cos(i1)]])
        m2 = numpy.array([[math.cos(theta)], [math.sin(theta)], [0]])   # 矩阵2
        xyz = r * numpy.dot(m1, m2)  # 根据公式矩阵点乘再数乘
        coord.append(xyz)
    return coord

site_3 = Read_N_file("D:\\Data\\卫星导航数据\\测点3\\20221201\\rinex\\3291471335I.22N")  # 读取测站3 N文件
site_12 = Read_N_file("D:\\Data\\卫星导航数据\\测点12\\测点12\\20221201\\rinex\\3291472335I.22N")  # 测站12 N文件
coordinate_3 = position(site_3)  # 站3观测卫星坐标计算结果
coordinate_12 = position(site_12)   # 站12观测卫星坐标计算结果

# 2、伪距单点定位
# 观测文件初始接收机坐标
A_xyz_0 = [-2615586.3441, 4732731.8405, 3371103.4148]  # 测站3模糊坐标
Ax0 = A_xyz_0[0]
Ay0 = A_xyz_0[1]
Az0 = A_xyz_0[2]
B_xyz_0 = [-2615466.3709, 4732800.5475, 3371105.5093]  # 测站12模糊坐标
Bx0 = B_xyz_0[0]
By0 = B_xyz_0[1]
Bz0 = B_xyz_0[2]
global x0   # 声明全局变量用于循环运算
global y0
global z0
site3_slite = ['G04', 'G16', 'G22', 'G26', 'G31']  # 5颗卫星，站3
site12_slite = ['G03', 'G04', 'G16', 'G22', 'G26', 'G27', 'G31']  # 7颗卫星，站7
C1c = [[22041120.015, 20400419.747, 21961188.555, 20504282.213, 21905730.086],
       [24141392.508, 21930913.665, 20290241.262, 21851108.433, 20394175.445, 21744151.797, 21795706.196]]  # 观测文件中伪距集合，[0]为站3，[1]为站12
site3_xyz_i = [coordinate_3[1], coordinate_3[2], coordinate_3[3], coordinate_3[4], coordinate_3[6]]  # 测站3卫星坐标其中G03卫星存在两个时刻：10：0：0和9:59:44，下标为0，暂不取用
# 测站12卫星坐标
site12_xyz_i = [coordinate_12[0], coordinate_12[1], coordinate_12[2], coordinate_12[3], coordinate_12[4],coordinate_12[5], coordinate_12[6]]

# 计算接收机坐标函数
def receiver_xyz(site_slite, n_site):
    count_number = 0  # 统计迭代次数
    c = 299792458  # 光速m/s
    while True:
        count_number += 1   # 迭代次数每轮加一
        list_A = []
        list_L = []
        global x0
        global y0
        global z0
        if count_number == 1:   # 防止初始坐标一直被赋值
            # 初始坐标赋值
            if len(site_slite) > 5:
                x0 = Bx0    # 测站12模糊坐标
                y0 = By0
                z0 = Bz0
            else:
                x0 = Ax0    # 测站3模糊坐标
                y0 = Ay0
                z0 = Az0
        for i in range(len(site_slite)):
            # 卫星坐标赋值
            if len(site_slite) > 5:  # 有效卫星数量大于5为测站12
                x_j = site12_xyz_i[i][0]
                y_j = site12_xyz_i[i][1]
                z_j = site12_xyz_i[i][2]
            else:   # 有效卫星数量5为测站3
                x_j = site3_xyz_i[i][0]
                y_j = site3_xyz_i[i][1]
                z_j = site3_xyz_i[i][2]
            distance0 = numpy.sqrt(numpy.power((x_j - x0), 2) + numpy.power((y_j - y0), 2) + numpy.power((z_j - z0), 2))  # 初始距离
            t = t1  # 要计算计算时间（转化为周内秒）
            toc = float(n_site[i]['t0'])    # 卫星参考时间
            tk = t - toc    # 归化时间
            delta_t = float(n_site[i]['卫星钟偏差']) + float(n_site[i]['卫星钟漂移']) * tk + float(n_site[i]['卫星钟漂移速度']) * math.pow(tk, 2)  # 卫星时钟误差
            lj = (x_j - x0) / distance0     # 伪距对X,Y,Z偏导l,m,n
            mj = (y_j - y0) / distance0
            nj = (z_j - z0) / distance0
            list_A_item = [float(lj), float(mj), float(nj), -c]     # A矩阵元素
            if len(site_slite) == 5:    # 区分测站3和测站12输入数据，测站3有效观测卫星5颗，测站12有效观测卫星7颗
                list_L_item = [distance0 - float(C1c[0][i]) - c * delta_t]  # 测站3观测文件伪距C1c
            else:
                list_L_item = [distance0 - float(C1c[1][i]) - c * delta_t]  # 测站12观测文件伪距C1c
            list_A.append(list_A_item)      # A矩阵元素加入用于构建矩阵的list中
            list_L.append(list_L_item[0])   # L矩阵元素加入用于构建矩阵的list中
        # 构建矩阵
        A = numpy.matrix(list_A)  # A
        L = numpy.matrix(list_L)    # L
        N = numpy.dot(A.T, A)   # A的转置乘A
        U = numpy.dot(A.T, L)   # A的转置乘L
        Target = numpy.dot(numpy.linalg.inv(N), U)  # 目标函数，一个1*4的矩阵存储deltaX,deltaY,deltaZ,deltaT
        x0 += Target[0, 0]  # 接收机坐标新赋值
        y0 += Target[1, 0]
        z0 += Target[2, 0]
        print("当前循环次数：", count_number)  # 统计循环次数
        print("当前接收机坐标：", "x：", x0, "y：", y0, "z：", z0)
        # 收敛条件，残差很小时候
        if (Target[0, 0] ** 2 + Target[1, 0] ** 2 + Target[2, 0] ** 2) < 0.000001 or count_number >= 50:
            A_XYZ = [Target[0, 0], Target[1, 0], Target[2, 0], Target[3, 0]]
            break
    return A_XYZ

# 主体部分
site3_result = receiver_xyz(site3_slite, site_3)  # 测站3信息
site12_result = receiver_xyz(site12_slite, site_12)     # 测站12位置
x3 = site3_result[0] + Ax0  # 测站3坐标
y3 = site3_result[1] + Ay0
z3 = site3_result[2] + Az0
x12 = site12_result[0] + Bx0    # 测站12坐标
y12 = site12_result[1] + By0
z12 = site12_result[2] + Bz0
Delta_t = [site3_result[3], site12_result[3]]   # 接收机钟差
length_true = 137.9468  # 以CGO2软件计算结果作为真值
print('测站3', 'x:', x3, 'y:', y3, 'z:', z3)
print('测站12', x12, y12, z12)
print('测站3接收机钟差:', Delta_t[0], '测站12接收机钟差：', Delta_t[1])
distance = math.dist((x3, y3, z3), (x12, y12, z12))     # 基线长为测站3与测站12欧式距离
print('基线解算：', distance, 'm')   # 基线结算结果
print('软件计算结果:', length_true, '\n差值为：', format(abs(distance - length_true), '.5f'), 'm')  # 小数保留5位有效数字
