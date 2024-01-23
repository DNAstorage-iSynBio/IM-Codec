# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 17:17:19 2021

@author: Wei Qiang
"""

from struct import pack
from md5hash import scan


met_tag = "Z"
nt_dic = {"A": 0, "C": 1, "G": 2, "T": 3}
nt_bin_dic = {
				4:[{"A":"11", "T":"10", "C":"01", "G":"00"}, 1],
				3: [{"AA":"000", "AT":"001", "AC":"010", "TT":"011", "TC":"101", "CC":"110", "CA":"111", "TA":"100"},2],
				2: [{"A":"1", "G":"0"},1],
				16: [{'A': '0000', 'B': '0001', 'C': '0010', 'D': '0011', 'E': '0100', 'F': '0101', 'G': '0110', 'H': '0111', 'I': '1000', 'J': '1001', 'K': '1010', 'L': '1011', 'M': '1100', 'N': '1101', 'O': '1110', 'P': '1111'}, 1]
				}


def removeArtificial(arti_info_seq_path):
	with open(arti_info_seq_path) as f1:
		arti_info_seq = f1.read().strip()

	reserve_num = arti_info_seq.count("Y")

	arti_decode_seq_path = arti_info_seq_path+".key"
	with open(arti_decode_seq_path) as f2:
		arti_decode_seq = f2.read().strip()
		reserve_seq, _tmp = arti_decode_seq.split("X")
		pos_seq, remove_seq = _tmp.split("Y")


	# Z index get
	reserve_pos_list = []
	for i in range(len(arti_info_seq)):
			if arti_info_seq[i] == "Y":
				reserve_pos_list.append(i)

	# remove_pos
	count_unit_len = int(len(reserve_seq)/reserve_num)
	remove_tag_num = int(len(remove_seq)/count_unit_len)
	remove_pos_list = []

	pos_unit_len = int(len(pos_seq)/remove_tag_num)
	for i in range(remove_tag_num):
		_pos_seq = pos_seq[pos_unit_len*i: pos_unit_len*(i+1)]
		pos = sijinzhiToShijinzhi(_pos_seq)
		reserve_pos_list.append(pos)

	# recover
	decode_seq_reverse_list = []
	arti_info_seq_list = list(arti_info_seq)

	tag_pos_list = reserve_pos_list + remove_pos_list
	tag_pos_list.sort(reverse=True)
	for pos in tag_pos_list:
		if arti_info_seq_list[pos] == "Y":
			arti_info_seq_list[pos] = "Z"
			decode_seq_reverse_list += [reserve_seq[-1*count_unit_len:]]
			reserve_seq = reserve_seq[:-2]

		else:
			arti_info_seq_list = arti_info_seq_list[:pos] + ["Z"] + arti_info_seq_list[pos:]
			decode_seq_reverse_list += [remove_seq[-1*count_unit_len:]]
			remove_seq = remove_seq[:-2]

	decode_seq = "".join(decode_seq_reverse_list[::-1])

	info_seq = "".join(arti_info_seq_list).replace("X", "ZZ")

	return info_seq, decode_seq










def transTag(seq:str, coding_nt_num:int) -> str:
	if coding_nt_num == 2:
		seq = seq.replace("C", "Z").replace("T","Z")
	if coding_nt_num == 3:
		seq = seq.replace("G","Z")

	return seq

def findEnd(seq: str) -> str:
	'''
	description: 将结束字符$确定位置并还原，结束字符由seq中 甲基化碱基连续个数第二多+1 表示，位置从第一多的第一个开始
	param seq: 信息序列
	return new_seq: 找到$符号并恢复的信息序列
	return decode_unit: 解码序列中，以几位为一个单元
	'''
	#甲基化碱基坐标
	met_list, len_list = [], []
	set_met_list = [0, 0]
	flag = False

	for pos in range(len(seq)):
		if seq[pos] == met_tag:
			if flag == False:
				flag = True
				set_met_list[0] = pos
		else:
			if flag == True:
				flag = False
				set_met_list[1] = pos
				len_list.append(set_met_list[1] - set_met_list[0])
				met_list.append(set_met_list)
				set_met_list = [0, 0]

	if flag == True:
		set_met_list[1] = pos + 1
		set_met_list = set_met_list[1] - set_met_list[0]
		met_list.append(set_met_list)
	len_uniq_list = list(set(len_list))
	len_uniq_list.sort()

	# 如果倒数第一位比倒数第二位多1，则$符号与序列中Z无连接，否则将最大数进行分解，满足新list中倒数第一位比倒数第二位多1
	if len_uniq_list[-1] != len_uniq_list[-2] + 1:

		for i in range(1, len_uniq_list[-1]):
			split_list = [i, len_uniq_list[-1]-i]
			len_uniq_list_tmp = len_uniq_list[:-1] + split_list
			# len_uniq_list_tmp = list(set(len_uniq_list_tmp))
			len_uniq_list_tmp.sort()

			if (len_uniq_list_tmp[-1] == len_uniq_list_tmp[-2] + 1) and (len_uniq_list_tmp[-1] in split_list):
				met_group_num = len(met_list)
				end_num = len_uniq_list_tmp[-1]
				break
	else:
		met_group_num = (len(met_list) - 1)
		end_num = len_uniq_list[-1]
	# print(len_uniq_list_tmp)

	# end_num = len_uniq_list[-1] - (len_uniq_list[-2] + 1)
	end_index = len_list.index(len_uniq_list[-1])
	end_list = met_list[end_index]

	new_seq = seq[:end_list[0]] + "$" + seq[end_list[0] + end_num:]

	return new_seq, met_group_num

def sijinzhiToShijinzhi(seq: str) -> int:
	'''
	description: 将碱基序列代表的四进制数转换为十进制
	'''
	sijinzhi_num_list = [nt_dic[i] for i in seq]
	shijinzhi_num = 0
	for i in range(len(sijinzhi_num_list)):
		shijinzhi_num += sijinzhi_num_list[i]*4**(len(sijinzhi_num_list)-i-1)



	return shijinzhi_num

def combineSeq(decode_seq, info_seq, met_group_num):
	# 得到解码序列所代表的十进制数字列表，即每个重复单元的重复次数
	decode_unit = int(len(decode_seq)/met_group_num)
	decode_list = []
	pos = 0
	while True:
		sijinzhi_seq = decode_seq[pos: pos + decode_unit]
		if sijinzhi_seq == "":
			break

		shijinzhi_num = sijinzhiToShijinzhi(sijinzhi_seq)
		decode_list.append(shijinzhi_num)
		pos += decode_unit
	combine_seq = ""
	decode_index, pos1 = 0, 0
	while True:
		try:
			info_seq[pos1]
		except:
			break
		if info_seq[pos1] == met_tag:
			for pos2 in range(pos1, len(info_seq)):
				if info_seq[pos2] != met_tag:
					break
			_seq_unit = info_seq[pos2: 2*pos2 - pos1]*decode_list[decode_index]
			combine_seq += _seq_unit
			pos1 = 2*pos2 - pos1
			decode_index += 1
		else:
			combine_seq += info_seq[pos1]
			pos1 += 1
	return combine_seq

def decodeBwt(last_seq):
	last_seq = list(last_seq)
	first_seq = [i for i in last_seq]
	first_seq.sort()


	recovery_seq = []
	seq_index = 0
	while True:
		_last_nt = last_seq[seq_index]
		recovery_seq.append(_last_nt)
		if _last_nt == "$":
			break

		# 相同碱基在first_seq和last_seq中的相对位置保持不变
		i = last_seq[:seq_index].count(_last_nt)
		seq_index = first_seq.index(_last_nt) + i
	recovery_seq.reverse()
	recovery_seq = "".join(recovery_seq[1:])
	return recovery_seq

def transBin(seq, coding_nt_num):
	bin_str = ""
	step = nt_bin_dic[coding_nt_num][1]
	for i in range(0, len(seq), step):
		_unit = seq[i: i+step]
		bin_str += nt_bin_dic[coding_nt_num][0][_unit]
	# bin_str = [nt_bin_dic[coding_nt_num][i] for i in seq]
	# bin_str = "".join(bin_str)

	## remove zero-padding
	zer_len = len(bin_str)%8
	if zer_len != 0:
		bin_str = bin_str[:-1*zer_len]

	return bin_str

# def transInfo(bin_str):
# 	bytes_list = []
# 	i = 0
# 	while True:
# 		if i >= len(bin_str):
# 			break

# 		_bin = int(bin_str[i: i+8], 2)
# 		bytes_list.append(_bin)
# 		i += 8
# 	info_str = bytes(bytes_list).decode("utf-8")
# 	# print("The information is: \n{0}".format(info_str))
# 	return info_str

def main(info_seq_ori, decode_seq, coding_nt_num,  output_path):
	info_seq_ori = transTag(info_seq_ori, coding_nt_num)
	info_seq, met_group_num = findEnd(info_seq_ori)
	combine_seq = combineSeq(decode_seq, info_seq, met_group_num)
	recovery_seq = decodeBwt(combine_seq)

	try:
		bin_str = transBin(recovery_seq, coding_nt_num)
		saveResult(bin_str, output_path)
		# info_str = transInfo(bin_str)
	except:
		info_str = "Error"
	# return info_str

def getSequence(info_path):
	key_path = info_path + ".key"
	 # info seq
	with open(info_path) as f:
		info_seq = f.read().strip()

	 # key seq
	with open(key_path) as f:
		key_seq = f.read().strip()

	return key_seq, info_seq

def saveResult(bin_str, result_path):
	with open(result_path, "wb") as f:
		i = 0
		while True:
			if i >= len(bin_str):
				break

			_bin = int(bin_str[i: i+8], 2)
			_a = pack('B', _bin)
			f.write(_a)
			i += 8


if __name__ == "__main__":

	import argparse
	def linxCommand():
		parser = argparse.ArgumentParser()
		parser.add_argument("i", type=str, help = "information path")
		parser.add_argument("o", type=str, help = "output path")
		parser.add_argument("-n", "--coding_nt_num", type=int, default=4, help = "The number of coding bases used, default=[4]")
		args = parser.parse_args()
		return args.i, args.o, args.coding_nt_num

	input_path, output_path, coding_nt_num = linxCommand()
	key_seq, info_seq = getSequence(input_path)
	info_str = main(info_seq, key_seq, coding_nt_num, output_path)
	# saveResult(info_str, output_path)










	# # 读文件获取序列
	# ## 人工操作部分需去除，如果有

	# arti_info_seq_path = "D:/BaiduSyncdisk/中科院先进院合成所/项目/2021_05_24-甲基化修饰DNA存储/03_测试数据/test/encode/JUNE9.txt.arti_e.seq"

	# info_seq, decode_seq = removeArtificial(arti_info_seq_path)

	## 读文件中的序列直接解码
	# info_seq_path = "D:/BaiduSyncdisk/中科院先进院合成所/项目/2021_05_24-甲基化修饰DNA存储/03_测试数据/test/encode/test.txt_e.seq"
	# decode_seq_path = info_seq_path +".key"

	# with open(info_seq_path) as f:
	# 	info_seq = f.read().strip()
	# with open(decode_seq_path) as f:
	# 	decode_seq = f.read().strip()

	# ## 直接写出序列
	# info_seq = "CTTTAATCCTTAAATTTAAAAATAAAAATATTCCCCAATTTATGGGGGTATAGCGATTCAGTATATTTAAACACTTACGGGGTAAATCA"
	# # # 真实 key
	# # decode_seq = "CGCTCTAT"

	# coding_nt_num = 4
	# output_path = "D:/BaiduSyncdisk/中科院先进院合成所/项目/2021_05_24-甲基化修饰DNA存储/03_测试数据/test/decode/test.txt"
	# info_str = main(info_seq, decode_seq, coding_nt_num, output_path)
	# print(scan(output_path))

	# # test key
	# import time
	# from met_dna_storage import siJinZhi

	# count = 4**10
	# _s = time.time()
	# for num_key in range(count):
	# 	decode_seq = siJinZhi(num_key)
	# 	# print(decode_seq)

	# 	try:
	# 		# 公共部分
	# 		coding_nt_num = 3
	# 		output_path = "D:/中科院先进院合成所/项目/2021_05_24-甲基化修饰DNA存储/10_加密/20230714/test.txt"
	# 		info_str = main(info_seq, decode_seq, coding_nt_num, output_path)
	# 	except:
	# 		pass
	# _e = time.time()
	# print("{} loop time spend:\t{}".format(count, _e - _s))

	# # # print(scan(output_path))

