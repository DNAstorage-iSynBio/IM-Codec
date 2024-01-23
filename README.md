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

```
python met_dna_storage.py input/test.txt output/text.dna -n 4 -m e > output/text.dna.log
```

## Help

Any advise for common problems or issues.

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

## Authors

Wei Qiang


## Version History

* Initial Release

## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Cite this program

*Xiaoluo Huang†*, Zhaohua Hou†, Wei Qiang†, Honglei Wang†, Xiangxiang Wang, Xiaoxu Chen, Xin Hu, Junbiao Dai*, Lingjun Li*, Guanghou Zhao*. **Towards next-generation DNA encryption via an expanded genetic system**. (Submitted to Journal)





