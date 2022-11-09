# -*- coding: utf-8 -*-
"""
Created on Tue May 25 11:17:33 2021

@author: Wei Qiang
"""

from PIL import Image 
import time, os

nt_dic = {"0": "A", "1": "C", "2": "G", "3": "T"}
met_tag = "Z"

def siJinZhi(x:int) -> str:
	'''
	description: 十进制数转为四进制，并根据对应关系转换为碱基序列
	param x: 十进制数
	return: 对应四进制数的碱基序列
	'''
	
	sijinzhi_list = []
	while x > 3:
		sijinzhi_list.append(str(x % 4))
		x //= 4
	sijinzhi_list.append(str(x))

	sijinzhi_list = [nt_dic[i] for i in sijinzhi_list]
	nt_seq = "".join(reversed(sijinzhi_list))
	return nt_seq

def trans_bin(string: str):
	'''
	description: 明文字符串转为01比特字符串

	param string: 明文字符串
	return: 01比特字符串
	'''
	bin_string = ""
	for _str in string:
		str_by = bytes(_str, encoding = "utf-8")
		# print(str_by)
		for str_b in str_by:
			str_bin = bin(str_b)
			# print(str_bin)
			str_bin_s = str(str_bin)[2:]
			str_bin_s = "0"*(8 - len(str_bin_s)) + str_bin_s
			bin_string += str_bin_s
	# print(len(bin_string))
	return bin_string

def transImageColor(image_path, image_type):
	'''
	description: 输入彩色图片，指定存储的类型【彩色，灰度，二值灰度】，程序中自动转换
	'''

	def intToBin(num):
		str_bin = bin(num)[2:]
		str_bin = "0"*(8-len(str_bin)) + str_bin
		return str_bin

	im=Image.open(image_path)
	width, height=im.size

	binary_str = ""
	for i in range(height):
		for j in range(width):
			color=im.getpixel((j,i))
			if image_type == "color":
				color_r, color_g, color_b = intToBin(color[0]), intToBin(color[1]), intToBin(color[2])
				binary_str += (color_r + color_g + color_b)
			elif image_type == "gray":
				color = int((color[0] + color[1] + color[2])/3)
				color = intToBin(color)
				binary_str += color
			elif image_type == "black_white":
				color = int((color[0] + color[1] + color[2])/3)
				if color < 128:
					binary_str += "0"
				else:
					binary_str += "1"

		# 添加换行符的二进制码
		binary_str += "00001010"

	return binary_str

def transImage(image_path, out_path):
	'''
	description: 将图片转换为二进制01字符串，仅黑白两色，一个像素点用一个0或1表示，其中0代表白色点，1代表黑色点

	param impage_path: 输入图片的路径
	param out_path: 可选参数，转为01字符串的结果保存为文本文件

	return binary_str: 二进制01字符串结果
	return (width, height): 图片的宽、高信息
	'''

	im=Image.open(image_path)
	width, height=im.size

	f = open(out_path, "w")
	binary_str = ""
	for i in range(height):
		for j in range(width):
			color=im.getpixel((j,i))
			colorsum=color[0]+color[1]+color[2]
			if colorsum < 256:
				binary_str += "0"
				f.write("0")
			else:
				binary_str += "1"
				f.write("1")
		f.write("\n")
	f.close()

	return binary_str, (width, height)

def transNtSeq(binary_str: str, coding_nt_num: int =4) -> str:
	'''
	description: 将二进制01字符串转换为碱基序列，对应关系为bin_nt_dic

	param binary_str: 二进制01字符串

	return nt_seq: 碱基序列
	'''
	if coding_nt_num == 4:
		bin_nt_dic = {"11":"A", "10":"T", "01":"C", "00":"G"}
		unit_len = 2
	elif coding_nt_num == 2:
		bin_nt_dic = {"1":"A", "0":"G"}
		unit_len = 1
	elif coding_nt_num == 3:
		bin_nt_dic = {"000":"AA", "001":"AT", "010":"AC", "011":"TT", "100":"TA", "101":"TC", "110":"CC", "111":"CA"}
		unit_len = 3

	n = 0
	nt_seq = ""
	while True:
		_unit = binary_str[n: n+unit_len]
		if _unit == "":
			break
		if len(_unit) < unit_len:
			print("Zero padding: {}\n".format(unit_len-len(_unit)))

			_unit += (unit_len-len(_unit))*"0"
 
		nt_seq += bin_nt_dic[_unit]
		n += unit_len

	return nt_seq


def suffixTree(string: str) -> list:
	'''
	description: 使用后缀树找出字符串中连续重复的子串（STR）

	param string: 字符串，此处为碱基序列

	return repeat_list: [[repeat_unit, start, end, length, repeat_count, unit_len], ...]，重复子字符串的相关信息，end为末尾pos+1，对应list取区间规则
	'''

	repeat_list = []

	for pos1 in range(len(string)): 
		for pos2 in range(pos1 + 1, len(string)):
			step_len = int(pos2 - pos1)

			#到达序列末尾，无法取到足够步长的序列，终止此层循环
			if len(string[pos2: pos2 + step_len]) < step_len:
				break

			if string[pos1: pos1 + step_len] == string[pos2: pos2 + step_len]:
				if repeat_list != [] and (int(len(string[pos1: pos1 + step_len])/repeat_list[-1][5])*repeat_list[-1][0]
				 == string[pos1: pos1 + step_len]) and (pos1 == repeat_list[-1][1]):
					continue

				count = 2
				while True:
					if string[pos2 + (count - 1)*step_len: pos2 + count*step_len] == string[pos1: pos1 + step_len]:
						count += 1
					else:
						break
				_repeat_list = [string[pos1: pos1 + step_len], pos1, pos1 + count*step_len, count*step_len, count, step_len]
				# if count == 2:
				# 	# print(_repeat_list)
				# 	continue

				#新找到的连续重复子串区间与现有最后一个子串区间重复时，若其长度大于上一个子串，则对上一子串进行替换
				if repeat_list != [] and (repeat_list[-1][1] <= _repeat_list[1] < repeat_list[-1][2]):
					if (_repeat_list[3] > repeat_list[-1][3]):
						repeat_list[-1] = _repeat_list
				else:
					repeat_list.append(_repeat_list)

	return repeat_list


def compressSeq(nt_seq: str):
	'''
	description: 识别序列中单一碱基的串联重复，并进行压缩

	param nt_seq: 碱基序列

	return comp_seq[1:]: 压缩后的序列
	return count_list[1:]: 序列中，碱基对应的计数信息
	'''

	comp_seq, count_list = "0", []
	n = 1

	for i in nt_seq:
		if i != comp_seq[-1]:
			count_list.append(n)
			comp_seq += i
			n = 1
		else:
			n += 1
	count_list.append(n)

	return comp_seq[1:], count_list[1:]


def outPutSeq(comp_seq, count_list):
	decode_seq_list, info_seq = [], ""
	set_len = 0
	for i in range(len(count_list)):
		_count = count_list[i]
		_nt = comp_seq[i]

		if _count == 1:
			info_seq += _nt
		else:
			_count_seq = siJinZhi(_count)
			decode_seq_list.append(_count_seq)
			if len(_count_seq) > set_len:
				set_len = len(_count_seq)

			info_seq += met_tag + _nt

	decode_seq_list = [nt_dic["0"]*(set_len - len(i)) + i for i in decode_seq_list]
	decode_seq = "".join(decode_seq_list)

	return decode_seq, info_seq, 


def bwt(nt_seq):
	'''
	description: BWT算法转换序列
	params nt_seq: 碱基序列
	return sort_list: 转换之后的矩阵
	'''
	nt_seq += "$"
	sort_list = [nt_seq]
	for i in range(len(nt_seq) - 1):
		nt_seq = nt_seq[-1] + nt_seq[:-1]
		#将新序列插入到有序list中，保持原有顺序

		pos1, pos2 = 0, len(sort_list) - 1
		while True:
			if pos2 == pos1 + 1:
				index1, index2 = pos1, pos2
			else:
				index1, index2 = round((pos1 + pos2)/2), round((pos1 + pos2)/2) + 1
			# print(nt_seq, sort_list, index1, index2, pos1, pos2) #bkp
			if index2 == len(sort_list):
				if nt_seq <= sort_list[index1]:
					sort_list = sort_list[:index1] + [nt_seq] + sort_list[index1:]
				else:
					sort_list.append(nt_seq)
				break
			else:
				if sort_list[index1] <= nt_seq and (nt_seq <= sort_list[index2]):
					sort_list = sort_list[:index2] + [nt_seq] + sort_list[index2:]
					break
				elif nt_seq < sort_list[index1]:
					pos2 = index1
				else:
					pos1 = index2
	last_seq = [i[-1] for i in sort_list]
	last_seq = "".join(last_seq)

	return sort_list, last_seq


def filterRepeatList(repeat_list):
	'''
	description: 过滤后缀树算法找到的串联重复序列，去掉压缩后总序列长度大于等于压缩前序列长的部分

	param repeat_list: 后缀树算法返回的串联重复序列的列表

	return repeat_list: 过滤之后的串联重复序列的列表
	return set_len: decode_seq中，以set_len长度为一个单元表示一组串联重复的重复个数
	return max_unit_len: 最长的重复单元的位数
	'''

	max_unit_len = max([i[5] for i in repeat_list])
	set_len = len(siJinZhi(max([i[4] for i in repeat_list])))

	# 连续两次过滤后的repeat_list相等时，跳出循环，返回参数
	while True:
		# print(max_unit_len, set_len, repeat_list,'\n') #bkp
		repeat_list_new = []
		for i in range(len(repeat_list)):
			repeat_unit, start, end, length, repeat_count, unit_len = repeat_list[i]
			if set_len + unit_len*2 < length:
				repeat_list_new.append(repeat_list[i])
		if repeat_list_new == repeat_list:
			break
		else:
			repeat_list = repeat_list_new
			max_unit_len = max([i[5] for i in repeat_list])
			set_len = len(siJinZhi(max([i[4] for i in repeat_list])))
			# print(max([i[4] for i in repeat_list])) #bkp
	return repeat_list, set_len, max_unit_len


# def filterRepeatList(repeat_list):

# 	[i.append(len(siJinZhi(i[4]))) for i in repeat_list]


# 	set_len_list = list(set([i[6] for i in repeat_list]))
# 	set_len_list.sort(reverse = True)
# 	decrease_nt_num_list = []
# 	for set_len in set_len_list:
# 		decrease_nt_num = 0
# 		for _repeat_list in repeat_list:
# 			repeat_count, unit_len, _set_len = _repeat_list[4], _repeat_list[5], _repeat_list[6]
# 			if _set_len > set_len:
# 				continue
# 			_devrease_nt_num = ((repeat_count-2)*unit_len- set_len)
# 			if _devrease_nt_num < 0:
# 				continue
# 			decrease_nt_num += _devrease_nt_num
# 		decrease_nt_num_list.append(decrease_nt_num)
# 	max_index = decrease_nt_num_list.index(max(decrease_nt_num_list))
# 	set_len = set_len_list[max_index]

# 	repeat_list_new = []
# 	for _repeat_list in repeat_list:
# 		repeat_count, unit_len, _set_len = _repeat_list[4], _repeat_list[5], _repeat_list[6]
# 		if _set_len > set_len:
# 			continue

# 		_devrease_nt_num = ((repeat_count-2)*unit_len- set_len)
# 		if _devrease_nt_num < 0:
# 			continue
# 		repeat_list_new.append(_repeat_list)

# 	max_unit_len = max([i[5] for i in repeat_list_new])

# 	return repeat_list_new, set_len, max_unit_len





def controlGC(nt_seq, win=50):
	'''
	for 2 baseUsed
	'''
	new_nt_seq = ""
	for i in range(0, len(nt_seq), win):
		_nt_seq = nt_seq[i: i+win]
		tag_num = _nt_seq.count(met_tag)

		g_num = _nt_seq.count("G")
		c_num = int(0.5*len(_nt_seq)-g_num)

		if c_num <= 0:
			_nt_seq = _nt_seq.replace(met_tag, "T")

		elif c_num < tag_num:
			_nt_seq = _nt_seq.replace(met_tag, "C", c_num)
			_nt_seq = _nt_seq.replace(met_tag, "T")

		else:
			_nt_seq = _nt_seq.replace(met_tag, "T")

		new_nt_seq += _nt_seq

	return new_nt_seq


def regularMain(nt_seq, coding_nt_num=4):
	

	# binary_str, size_set = transImage(image_path, out_path)
	# nt_seq = transNtSeq(binary_str)

	#bwt
	print("Start: the step of BWT")
	s = time.time()
	sort_seq_matrix, trans_seq = bwt(nt_seq)
	e = time.time()
	print("Done! Time: {0} min.\n".format((e - s)/60))
	# trans_seq = None

	#后缀树
	print("Start: the step of suffix Tree")
	s = time.time()
	repeat_list_ori = suffixTree(trans_seq)
	e = time.time()
	print("Done! Time: {0} min.\n".format((e - s)/60))

	#过滤无法有压缩效果的串联重复序列（压缩后的序列总长大于等于压缩前的序列）
	print("Start: the step of filter repeat list")
	s = time.time()
	repeat_list, set_len, max_unit_len = filterRepeatList(repeat_list_ori)
	e = time.time()
	print("Done! Time: {0} min.\n".format((e - s)/60))

	decode_seq_list, info_seq = [], ""
	# set_len = 0
	# max_unit_len = 0

	## 得到编码序列和信息序列
	break_point_index = 0
	for _repeat_list in repeat_list:
		repeat_unit, start, end, length, repeat_count, unit_len = _repeat_list[:6]
		# if unit_len > max_unit_len:
		# 	max_unit_len = unit_len
		_decode_seq = siJinZhi(repeat_count)
		decode_seq_list.append(_decode_seq)
		# if len(_decode_seq) > set_len:
		# 	set_len = len(_decode_seq)

		_info_seq = trans_seq[break_point_index: start] + unit_len*met_tag + trans_seq[start: start + unit_len]
		info_seq += _info_seq
		break_point_index = end
	info_seq += trans_seq[end:]

	# 编码序列
	print("The length of decode-sequence-unit is: {0}. ".format(set_len))
	decode_seq_list = [nt_dic["0"]*(set_len - len(i)) + i for i in decode_seq_list]
	decode_seq = "".join(decode_seq_list)

	#替换$符号为(set_len+1)个甲基化碱基
	tag_index = info_seq.index("$")
	info_seq = info_seq[:tag_index] + met_tag*(max_unit_len + 1) + info_seq[tag_index + 1:]
	if coding_nt_num==2:
		info_seq = controlGC(info_seq)
	elif coding_nt_num==3:
		info_seq = info_seq.replace(met_tag, "G")

	#压缩效率
	compressibility = round((len(decode_seq) + len(info_seq))/len(nt_seq)*100, 2)

	print("The compressibility is: {0}%.".format(compressibility))
	return decode_seq, info_seq, nt_seq, trans_seq, repeat_list, sort_seq_matrix, repeat_list_ori

def singleBaseCompressMain(nt_seq):
	'''
	description: 压缩序列中一个碱基为单元的串联重复
	param nt_seq: 碱基序列
	return:
	'''

	comp_seq, count_list = compressSeq(nt_seq)
	decode_seq, info_seq = outPutSeq(comp_seq, count_list)
	compressibility = round((len(decode_seq) + len(info_seq))/len(nt_seq)*100, 2)
	print("The compressibility is: {0}%.".format(compressibility))
	return decode_seq, info_seq

def transBin(input_path):
	binary_str = ""
	with open(input_path, "rb") as f:
		for line in f:
			for i in line:
				_bin_str = bin(i)[2:]
				_bin_str = "0"*(8-len(_bin_str))+_bin_str
				binary_str += _bin_str
	return binary_str

def readFile(input_path:str) -> str:
	with open(input_path) as f:
		file_str = f.read().strip()

	return file_str

def saveResult(out_path, decode_seq, info_seq):
	decode_path = out_path + ".key"
	with open(out_path, "w") as f:
		f.write(info_seq+"\n")

	with open(decode_path, "w") as f:
		f.write(decode_seq+"\n")

# def saveTmp(path, data):
# 	f = open(path, "w")
# 	f.write(data)
# 	f.close()
if __name__ == "__main__":
	# # image_path = "D:/中科院先进院合成所/项目/2021_05_24-甲基化修饰DNA存储/测试数据/向日葵250318.jpg"
	# # out_path = "D:/中科院先进院合成所/项目/2021_05_24-甲基化修饰DNA存储/output/太极.txt"
	# # binary_str, size_set = transImage(image_path, out_path)
	# # nt_seq = transNtSeq(binary_str)
	# # decode_seq, info_seq = singleBaseCompressMain(nt_seq)


	# # poetry = "Hello biology!"
	# poetry = "Hello world!"
	# # poetry = '''The cloud stood hunbly in a corner of the sky, the morning crowned it with splend our. '''
	# # poetry = "Information density is a further consideration while assessing this technique. Figure files have the lowest information density for each approach."


	# # input_path = "D:/中科院先进院合成所/项目/2021_05_24-甲基化修饰DNA存储/文字测试数据/input/背影-朱自清.txt"
	# # with open(input_path, encoding = "utf-8") as f:
	# # 	poetry = [i for i in f]
	# # poetry = "".join(poetry)
	
	# binary_str = trans_bin(poetry)

	# # input_path = "D:/中科院先进院合成所/项目/2021_05_24-甲基化修饰DNA存储/03_测试数据/加密相关测试数据/input/1.jpg"
	# # input_path = "D:/中科院先进院合成所/项目/2021_05_24-甲基化修饰DNA存储/11_图片/test.txt"
	# # input_path = "D:/中科院先进院合成所/项目/2021_05_24-甲基化修饰DNA存储/07_实验/input.txt"
	# # input_path = "D:/中科院先进院合成所/项目/2021_05_24-甲基化修饰DNA存储/12_bitcoin/praviteKey.txt"
	# # binary_str = transBin(input_path)
	# # binary_str = transImageColor(image_path, "black_white")
	# coding_nt_num=2
	# nt_seq = transNtSeq(binary_str, coding_nt_num=coding_nt_num)
	# decode_seq, info_seq, nt_seq, trans_seq, repeat_list, sort_seq_matrix, repeat_list_ori = regularMain(nt_seq, coding_nt_num=coding_nt_num)
	# # saveResult(out_path, decode_seq, info_seq)

	###########################################################################################################################
	###########################################################################################################################
	import argparse
	def linxCommand():
		parser = argparse.ArgumentParser()
		parser.add_argument("i", type=str, help = "input path")
		parser.add_argument("o", type=str, help = "output path")
		parser.add_argument("-n", "--coding_nt_num", type=int, default=4, help = "The number of coding bases used, default=[4]")
		parser.add_argument("-bin", "--binary", type=bool, default=False, help="The input file is binary string or not, default=[False]. Bool type.")
		args = parser.parse_args()
		return args.i, args.o, args.coding_nt_num, args.binary

	input_path, out_path, coding_nt_num, is_bin = linxCommand()

	if is_bin == False:
		binary_str = transBin(input_path)
	else:
		binary_str = readFile(input_path)

	# print(type(binary_str)) #bkp
	nt_seq = transNtSeq(binary_str, coding_nt_num=coding_nt_num)
	decode_seq, info_seq, nt_seq, trans_seq, repeat_list, sort_seq_matrix, repeat_list_ori = regularMain(nt_seq)
	saveResult(out_path, decode_seq, info_seq)





	






