# -*- coding: utf-8 -*-
"""
Created on Tue May 25 11:17:33 2021

@author: Wei Qiang
"""

from PIL import Image 
import time, os
from md5hash import scan

nt_dic = {"0": "A", "1": "C", "2": "G", "3": "T"}
met_tag = "Z"

def siJinZhi(x:int) -> str:
	'''
	description: Decimal numbers are converted to quad and converted to base sequences according to correspondence.
	param x: Decimal number
	return: base_sequence
	'''
	
	sijinzhi_list = []
	while x > 3:
		sijinzhi_list.append(str(x % 4))
		x //= 4
	sijinzhi_list.append(str(x))

	sijinzhi_list = [nt_dic[i] for i in sijinzhi_list]
	nt_seq = "".join(reversed(sijinzhi_list))
	return nt_seq




def transNtSeq(binary_str: str, coding_nt_num: int =4) -> str:
	'''
	description: trasform binary tring into base sequence, mapping relationship: bin_nt_dic

	param binary_str: binary string

	return nt_seq: base sequence
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
	elif coding_nt_num == 8:
		bin_nt_dic = {"000":"A", "001":"B", "010":"C", "011":"D", "100":"E", "101":"F", "110":"G", "111":"H"}
		unit_len = 3
	elif coding_nt_num == 16:
		bin_nt_dic = {"0000": "A", "0001": "B", "0010": "C", "0011": "D", "0100": "E", "0101": "F", "0110": "G",
					  "0111": "H", "1000": "I", "1001": "J", "1010": "K", "1011": "L", "1100": "M", "1101": "N",
					  "1110": "O", "1111": "P"}
		unit_len = 4
		
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
	description: find the tandem repeat sub-sequence

	param string: base sequence

	return repeat_list: [[repeat_unit, start, end, length, repeat_count, unit_len], ...]，information of tandem repeat sub-sequence, end = last_pos+1
	'''

	repeat_list = []

	for pos1 in range(len(string)): 
		for pos2 in range(pos1 + 1, len(string)):
			step_len = int(pos2 - pos1)

			# Arriving at the end of a sequence, it is impossible to take a sequence with sufficient step length, terminate this layer loop
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

				'''
				When the newly found continuous repeat substring interval repeats with the existing last substring interval, 
				if its length is greater than the previous substring, replace the previous substring
				'''
				if repeat_list != [] and (repeat_list[-1][1] <= _repeat_list[1] < repeat_list[-1][2]):
					if (_repeat_list[3] > repeat_list[-1][3]):
						repeat_list[-1] = _repeat_list
				else:
					repeat_list.append(_repeat_list)

	return repeat_list

def suffixTreeEcryption(string: str) -> list:
	# record the single base repeat unit

	# add marker at the end 
	string += "?"
	s = 0
	repeat_list = []

	while True:
		n = 1
		if s == len(string)-1:
			break

		while True:
			base_s, base_e = string[s], string[s+n]
			if base_e != base_s:
				_repeat_list = [base_s, s, s+n, n, n, 1]
				repeat_list.append(_repeat_list)
				break
			else:
				n+= 1

		s += n

	return repeat_list







def bwt(nt_seq):
	'''
	description: BWT algorithm tranform the base sequence
	params nt_seq: base sequence
	return sort_list: the matirx after tramsformation
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
	description: filter the tandem sequencesm. remove the sequence that the length after compression is bigger than original.

	param repeat_list: tandem repeat sequences list

	return repeat_list: tandem repeat sequences list after filter
	return set_len: in decode_seq，the length of the number of unit repeating
	return max_unit_len
	'''

	max_unit_len = max([i[5] for i in repeat_list])
	set_len = len(siJinZhi(max([i[4] for i in repeat_list])))

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


def encryptFilterList(repeat_list, run = 4):
	# filter in encryption mode. choose the single base repeating sub-sequence
	# run: homopolymer can be set, default value=4 

	repeat_list_new = []
	for _repeat_list in repeat_list:
		if _repeat_list[-1] == 1:
			if _repeat_list[4] > run:
				repeat_list_new.append(_repeat_list)

	max_unit_len = 1
	set_len = len(siJinZhi(max([i[4] for i in repeat_list_new])))

	return repeat_list_new, set_len, max_unit_len


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



	return repeat_list

def regularMain(nt_seq, coding_nt_num=4, mode="c"):
	


	#bwt
	print("Start: the step of BWT")
	s = time.time()
	sort_seq_matrix, trans_seq = bwt(nt_seq)
	e = time.time()
	print("Done! Time: {0} min.\n".format((e - s)/60))
	# trans_seq = None

	print("Start: the step of suffix Tree")
	s = time.time()

	if mode == "c":
		repeat_list_ori = suffixTree(trans_seq)
	elif mode == "e":
		repeat_list_ori = suffixTreeEcryption(trans_seq)

	e = time.time()
	print("Done! Time: {0} min.\n".format((e - s)/60))

	print("Start: the step of filter repeat list")
	s = time.time()

	if mode == "c":
		repeat_list, set_len, max_unit_len = filterRepeatList(repeat_list_ori)
	elif mode == "e":
		repeat_list, set_len, max_unit_len = encryptFilterList(repeat_list_ori)
	e = time.time()

	print("Done! Time: {0} min.\n".format((e - s)/60))

	decode_seq_list, info_seq = [], ""
	# set_len = 0
	# max_unit_len = 0

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

	print("The length of decode-sequence-unit is: {0}. ".format(set_len))
	decode_seq_list = [nt_dic["0"]*(set_len - len(i)) + i for i in decode_seq_list]
	decode_seq = "".join(decode_seq_list)

	# change the "$" to (set_len+1)*met_tag
	tag_index = info_seq.index("$")
	info_seq = info_seq[:tag_index] + met_tag*(max_unit_len + 1) + info_seq[tag_index + 1:]
	if coding_nt_num==2:
		info_seq = controlGC(info_seq)
	elif coding_nt_num==3:
		info_seq = info_seq.replace(met_tag, "G")

	# compressibililty
	compressibility = round((len(decode_seq) + len(info_seq))/len(nt_seq)*100, 2)

	print("The compressibility is: {0}%.".format(compressibility))
	return decode_seq, info_seq, nt_seq, trans_seq, repeat_list, sort_seq_matrix, repeat_list_ori


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
	out_path += "_{}.seq".format(mode)
	decode_path = out_path + ".key"
	with open(out_path, "w") as f:
		f.write(info_seq+"\n")

	with open(decode_path, "w") as f:
		f.write(decode_seq+"\n")

# remove_index_list 指的是移除第几个Z, 不是在序列中的index
def artificialChoosePos(remove_index_list:list, info_seq:str, decode_seq: str):
	
	reserve_seq, pos_seq, remove_seq = "", "", ""

	arti_info_seq = info_seq.replace(met_tag*2, "X")
	tag_num = len(repeat_list)
	unit_len = int(len(decode_seq)/tag_num)


	pos_len = len(siJinZhi(len(arti_info_seq)))
	for i in range(tag_num):
		if i not in remove_index_list:
			arti_info_seq = arti_info_seq.replace("Z", "Y", 1)
			reserve_seq += decode_seq[unit_len*i: unit_len*(i+1)]
		else:
	
			tag_index = arti_info_seq.index("Z")
			arti_info_seq = arti_info_seq.replace("Z", "", 1)

			_pos_seq = siJinZhi(tag_index)
			_pos_seq = nt_dic["0"]*(pos_len-len(_pos_seq))+_pos_seq	
			pos_seq += _pos_seq

			remove_seq += decode_seq[unit_len*i: unit_len*(i+1)]

	arti_decode_seq = "{}X{}Y{}".format(reserve_seq, pos_seq, remove_seq)

	return arti_decode_seq, arti_info_seq






if __name__ == "__main__":

	import argparse
	def linxCommand():
		parser = argparse.ArgumentParser()
		parser.add_argument("i", type=str, help = "input path")
		parser.add_argument("o", type=str, help = "output path")
		parser.add_argument("-n", "--coding_nt_num", type=int, default=4, help = "The number of coding bases used, default=[4]")
		parser.add_argument("-m", "--mode", type=str, default="e", help="[e]ncryption, [c]ompression.")
		parser.add_argument("-bin", "--binary", type=bool, default=False, help="The input file is binary string or not, default=[False]. Bool type.")
		args = parser.parse_args()
		return args.i, args.o, args.coding_nt_num, args.binary, args.mode

	input_path, out_path, coding_nt_num, is_bin, mode = linxCommand()

	if is_bin == False:
		binary_str = transBin(input_path)
	else:
		binary_str = readFile(input_path)

	
	nt_seq = transNtSeq(binary_str, coding_nt_num=coding_nt_num)
	decode_seq, info_seq, nt_seq, trans_seq, repeat_list, sort_seq_matrix, repeat_list_ori = regularMain(nt_seq, coding_nt_num=coding_nt_num, mode=mode)
	saveResult(out_path, decode_seq, info_seq)

	file_size = os.path.getsize(input_path)*8
	print("{} density: {} bits/nt\n\n".format(os.path.basename(input_path), file_size/(len(info_seq)+len(decode_seq))))







	






