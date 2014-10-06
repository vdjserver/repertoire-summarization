#!/usr/bin/env python

def create_read_id_count_map(csv_file_path):
	reader=open(csv_file_path)
	effective_line_num=0
	read_count_map=dict()
	for line in reader:
		temp_line=line.strip()
		if(temp_line[0]=="#"):
			#ignore lines starting with #
			pass
		else:
			if(effective_line_num==0):
				#ignore header 
				pass
			else:
				pieces=temp_line.split('\t')
				read_id=pieces[0]
				count_str=pieces[1]
				count_val=int(count_str)
				read_count_map[read_id]=count_val
	reader.close()
	return read_count_map
		











