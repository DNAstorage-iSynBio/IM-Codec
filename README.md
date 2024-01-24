# IM-Codec DNA storage algorithm

IM-Codec algorithm encrypts data into a “key sequence” (KS) and an “information sequence” (IS).

## Description

The IM-codec algorithm process the running length in the encoded sequence to encrypt sequence information and improve storage density.

## Getting Started

### Dependencies



* Windows 10 or Linux
* Python 3.10


### Installing

* https://github.com/DNAstorage-iSynBio/IM-Codec.git

### Executing program demo

Encoding

```
python met_dna_storage.py input/test.txt encode/test -n 4 -m e > encode/test.log
```

Decoding

```
python met_dna_decode.py encode/test.seq decode/test.seq.txt -n 4 > decode/test.seq.log
```

## Help

Encoding

```
python met_dna_storage.py -h
```

```
usage: met_dna_storage.py [-h] [-n CODING_NT_NUM] [-m MODE] [-bin BINARY] i o

positional arguments:
  i                     input path
  o                     output path

optional arguments:
  -h, --help            show this help message and exit
  -n CODING_NT_NUM, --coding_nt_num CODING_NT_NUM
                        The number of coding bases used, default=[4]
  -m MODE, --mode MODE  [e]ncryption, [c]ompression.
  -bin BINARY, --binary BINARY
                        The input file is binary string or not, default=[False]. Bool type.
```

Decoding

```
python met_dna_decode.py  -h
```

```
usage: met_dna_decode.py [-h] [-n CODING_NT_NUM] i o

positional arguments:
  i                     information path
  o                     output path

optional arguments:
  -h, --help            show this help message and exit
  -n CODING_NT_NUM, --coding_nt_num CODING_NT_NUM
                        The number of coding bases used, default=[4]
```

## Code author

Wei Qiang


## Version History

* Initial Release

## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Cite this program

Xiaoluo Huang†*, Zhaohua Hou†, Wei Qiang†, Honglei Wang†, Xiangxiang Wang, Xiaoxu Chen, Xin Hu, Junbiao Dai*, Lingjun Li*, Guanghou Zhao*. **Towards next-generation DNA encryption via an expanded genetic system**. (Submitted to Journal)





