#!/bin/bash
make -f makefileInit
./Init
echo -e "\n"
echo $?
echo "---------------------------------------|"
echo "                                       |"
echo "                                       |"
echo "    FILE INITIALIZED                   |"
echo "                                       |"
echo "                                       |"
echo "---------------------------------------|"
make
./Main


