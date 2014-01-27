#!/usr/bin/python
import urllib2
from bs4 import BeautifulSoup
from collections import defaultdict
import pprint

def tree(): return defaultdict(tree)
def dicts(t): return {k: dicts(t[k]) for k in t}

my_hier=tree()
url="http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=genetable&species=human&group=IGHV"
print "USING HTML RETRIEVED FROM ",url
f = urllib2.urlopen(url)
#print f.read()
html=f.read()
#fh=open("index.php?section=LocusGenes&repertoire=genetable&species=human&group=IGHD",'r')
#html=""
#for line in fh:
#	html+=line
#fh.close()

soup=BeautifulSoup(html)
#print(soup.prettify())
table_div=soup.find(id="genetable")
#print "GENE TABLE DIV:"
#print table_div.prettify()
gene_table=table_div.find('table')
#print "GENE TABLE:"
#print gene_table
#gene_table=Beautiful
rows = gene_table.findAll('tr')
subgroup=""
gene_name=""
allele_name=""
group="IGHD"
for tr in rows:
	#print "\tEXPLORING A TR...."
	cols = tr.findAll('td')
	col_num=0
	for td in cols:
		#print "\t\tEXPLORING A TD (col_num=",col_num,")... type=",type(td)
		if('class' in td.attrs):
			#print "class in td.attrs"
			class_val=td.attrs['class']
			#print "The class value is ",class_val
			class_val_first=str(class_val[0])
			#print "class_val_first is ",class_val_first
			class_val_first_str=str(class_val_first)
			if(class_val_first_str=="subgroup_middle_note"):
				text = td.find(text=True)# + ';'
				if text is None:
					text="NONE;"
				else:
					text+=";"
				#print "GOT A SUBGROUP=",text
				subgroup=str(text)
			elif(class_val_first_str=="gene_note"):
				text = td.find(text=True)# + ';'
				if text is None:
					text="NONE;"
				else:
					text+=";"
				#print "GOT A GENE=",text
				gene_name=str(text)
			elif(class_val_first_str=="allele_note"):
				text = td.find(text=True)# + ';'
				if text is None:
					text="NONE;"
				else:
					text+=";"
				#print "GOT AN ALLELE=",text
				allele_name=str(text)
				my_hier[group][subgroup][gene_name][allele_name]


pp = pprint.PrettyPrinter()
pp.pprint(dicts(my_hier))
