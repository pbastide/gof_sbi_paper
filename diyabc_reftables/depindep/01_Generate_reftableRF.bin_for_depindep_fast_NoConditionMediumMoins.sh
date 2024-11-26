# Copyright (C) {2024} {GLM, PB, AE, JMM}
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# Bash script
# Right-click to open a MobaXterm terminal (virtual Linux environment for Windows)
# Make the script executable by running the command: chmod +x 01_Generate_reftableRF.bin_for_depindep_fast_NoConditionMediumMoins.sh
# Run the script with: ./01_Generate_reftableRF.bin_for_depindep_fast_NoConditionMediumMoins.sh
#
# REQUIREMENTS:
# - DIYABC Version 1.1.51 or above.
#   Please download from https://github.com/diyabc/diyabc/releases/tag/v1.1.51
#   and change variable "diyabc" below with the path to your executable file.
# - data file DSIM_indep_dep_4pop_SNPind_10indsPerPop.snp
# - header file headerRF_for_depindep_fast_NoConditionMediumMoins.txt (which has to be renamed headerRF.txt)
#
#!/bin/bash
diyabc=diyabc-RF-windows-v1.1.51.exe # path to diyabc executable
./${diyabc} -p ./ -n "t:32;c:1;s:1;f:f"
./${diyabc} -p ./ -R "ALL" -r 220000 -g 1000 -m -t 32
