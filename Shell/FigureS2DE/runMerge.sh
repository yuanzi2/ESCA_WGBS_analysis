#!/bin/sh
#$ -l h_vmem=80G

cd Data/Model2D3D_PMDs/PMDs
bedtools multiinter -i EAC_1.Model2D3D.rmBlackList.bed EAC_2.Model2D3D.rmBlackList.bed EAC_3.Model2D3D.rmBlackList.bed EAC_4.Model2D3D.rmBlackList.bed \
                       EAC_6.Model2D3D.rmBlackList.bed GEJ_1.Model2D3D.rmBlackList.bed GEJ_2.Model2D3D.rmBlackList.bed GEJ_3.Model2D3D.rmBlackList.bed \
                       GEJ_4.Model2D3D.rmBlackList.bed GEJ_5.Model2D3D.rmBlackList.bed GEJ_6.Model2D3D.rmBlackList.bed GEJ_7.Model2D3D.rmBlackList.bed \
                       -g ../../meta/hg38_length.bed -header  > ../EAC_PMDs_multiinter.bed

bedtools multiinter -i ESCC_1.Model2D3D.rmBlackList.bed ESCC_2.Model2D3D.rmBlackList.bed ESCC_3.Model2D3D.rmBlackList.bed ESCC_4.Model2D3D.rmBlackList.bed \
                       ESCC_5.Model2D3D.rmBlackList.bed ESCC_6.Model2D3D.rmBlackList.bed ESCC_7.Model2D3D.rmBlackList.bed ESCC_8.Model2D3D.rmBlackList.bed \
                       ESCC_9.Model2D3D.rmBlackList.bed ESCC_10.Model2D3D.rmBlackList.bed ESCC_11.Model2D3D.rmBlackList.bed ESCC_12.Model2D3D.rmBlackList.bed \
                       ESCC_13.Model2D3D.rmBlackList.bed ESCC_14.Model2D3D.rmBlackList.bed ESCC_15.Model2D3D.rmBlackList.bed ESCC_16.Model2D3D.rmBlackList.bed \
                       ESCC_17.Model2D3D.rmBlackList.bed ESCC_19.Model2D3D.rmBlackList.bed ESCC_20.Model2D3D.rmBlackList.bed ESCC_21.Model2D3D.rmBlackList.bed \
                       ESCC_22.Model2D3D.rmBlackList.bed  -g ../../meta/hg38_length.bed -header > ../ESCC_PMDs_multiinter.bed
