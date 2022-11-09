# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 17:17:19 2021

@author: Wei Qiang
"""

met_tag = "Z"
nt_dic = {"A": 0, "C": 1, "G": 2, "T": 3}
nt_bin_dic = {
				4:[{"A":"11", "T":"10", "C":"01", "G":"00"}, 1],
				3: [{"AA":"000", "AT":"001", "AC":"010", "TT":"011", "TC":"101", "CC":"110", "CA":"111", "TA":"100"},2],
				2: [{"A":"1", "G":"0"},1]
				}
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
	return bin_str

def transInfo(bin_str):
	bytes_list = []
	i = 0
	while True:
		if i >= len(bin_str):
			break

		_bin = int(bin_str[i: i+8], 2)
		bytes_list.append(_bin)
		i += 8
	info_str = bytes(bytes_list).decode("utf-8")
	# print("The information is: \n{0}".format(info_str))
	return info_str

def main(info_seq_ori, decode_seq, coding_nt_num):
	info_seq_ori = transTag(info_seq_ori, coding_nt_num)
	info_seq, met_group_num = findEnd(info_seq_ori)
	combine_seq = combineSeq(decode_seq, info_seq, met_group_num)
	recovery_seq = decodeBwt(combine_seq)
	
	try:
		bin_str = transBin(recovery_seq, coding_nt_num)
		info_str = transInfo(bin_str)
	except:
		info_str = "Error"
	# print(info_str)
	return info_str

def getSequence(info_path):
	key_path = args.i + ".key"
	 # info seq
	 with open(info_path) as f:
	 	info_seq = f.read().strip()

	 # key seq
	 with open(key_path) as f:
	 	key_seq = f.read().strip()

	 return key_path, info_seq

def saveResult(info_str, result_path):
	with open(result_path, "w") as f:
		f.write(info_str+"\n")

if __name__ == "__main__":

	import argparse
	def linxCommand():
		parser = argparse.ArgumentParser()
		parser.add_argument("i", type=str, help = "information path")
		parser.add_argument("o", type=str, help = "output path")
		parser.add_argument("-n", "--coding_nt_num", type=int, default=4, help = "The number of coding bases used, default=[4]")
		args = parser.parse_args()
		return args.i, args.o, args.coding_nt_num

		key_seq, info_seq = getSequence(args.i)
		info_str = main(info_seq, key_seq, args.coding_nt_num)
		saveResult(info_str, args.o)












	# info_seq_ori = "CTTTTACTAAGCTAAACCCATTCAATGCAAAACCCACCACTTGGGGCCCACTTCATATTAGGATTTATTTCTTTAAATCCCATTTATTTTTTAACCAATAATTCTTTTTGGGTTACTTACCAAAACTATCTAAACGACATACCCCTCCACCAACCCACCCTCCAAATGGGCCAATGACCCTTATACTACAAACTTTTCTTTTCGATGCATCCTTGCTTCCTCCTTATATCTCCTCTTCTACCAAAATATTCTCTTTTTATTTACTTAATAAGTGCTCTTCCCTCAAACCTTCTTTCTAAACTCCCAAATTTAACGACTCATATAATTTTATCACATATTTAAACTAAAATATTAATTACTTTATAACTCCTTTCTCATTTTACAA"
	# decode_seq = "CCCCCAATCTATCCCCCTCCGCCCGA"
	# coding_nt_num = 3
	# info_str = main(info_seq_ori, decode_seq, coding_nt_num)
