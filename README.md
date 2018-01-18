# clas12_analysis

Quick test to show hipo2root is working and output a few histograms to a root file.

Prerequisites:
* [cern root](https://root.cern.ch/)
* pyroot

### cpp
A c++ example can be found in the cpp folder.


To build:

```bash
make clean && make
```

To run:
```bash
./test /path/to/root/input.root /path/to/output.root
```

To run over multiple files:
```bash
./test '/path/to/root/input/*.root' /path/to/output.root
```

### python
A python example can be found in the python folder.

To run:
```
./main.py '/path/to/root/input/*.root'
```

### Java

Java is dumb, don't use it. If you must use java, there's no easy way, so copy and paste the example.java file into an eclipse environment with the working [coatjava](https://github.com/JeffersonLab/clas12-offline-software.git) libraries. I think I'm using version 5.0 (Sorry I mean 5a.0.x because there's an a? Version numbers usually have letters right?) but the naming scheme and version numbers are so messed up I don't know what any of it means.

If it doesn't work don't complain to me..... :information_desk_person:
