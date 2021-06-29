#!/bin/bash

#run=5225
#run=5460
run=5387
#run=5303
#run=5308

class_name="ProtonMassProduction_run"
class_namex=$class_name"X"

selector_name=$class_name$run
ana_name="makeproton_ana_"$class_name$run

header_name=$class_name$run".h"
tmp_header_name=$class_name${run}"_tmp.h"

class_code=$class_name${run}".C"
tmp_class_code=$class_name$run"_tmp.C"



echo $class_namex

#[1]Generate the file list for analysis

#[2]Generate the ana module [to get the dat structure of selected trees]
g++ makeproton_ana.cc `root-config --libs --cflags` -o $ana_name


#[3]Run the ana module (input can be changable if needed but still need compile to loop over the selected files)
./$ana_name list_data_prod4_run$run'.txt' $selector_name

#[4]Fix bugs in the generated makeclass module
sed '/Init(tree)\;/i if (tree-\>InheritsFrom(\"TChain\")) ((TChain\*)tree)-\>LoadTree(0);' $header_name > $tmp_header_name
mv $tmp_header_name $header_name

#[4.2]Fix bug of GetEntry function (not necessary process if we a template already, which is the case for the moment)
#sed 's/GetEntriesFast/GetEntries/g' ProtonSelector_run${run}.C > ProtonSelector_0_run${run}.C

#[4.3]copy an existing code and replace the string to the selected run
cp -prv $class_namex".C"  $tmp_class_code
sed 's/5387/'${run}'/g' $tmp_class_code > $class_code
rm -f $tmp_class_code

root_exe_str="root -b -q 'RunAna.C(\""$selector_name\"")'"

echo $root_exe_str" ......"
eval $root_exe_str

