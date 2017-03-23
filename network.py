# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 10:24:54 2017

@author: lu
"""

import networkx as nx
import pickle as p
import pylab as pl 
import matplotlib.pyplot as  plt

# id-namespaces
def id_namespaces():  
    with open('goid_namespaces_dic.txt','w') as f1:
        with open(r'G:\files\go.obo', 'r') as f:
            dic = {}
            for line in f:
			line_data=line.strip().split(':',1)
				#line_data=line.strip().split(':')
			if line_data[0]=='id':
				son_node=line_data[1].strip()
					#son_node=line_data[1]
			if line_data[0]=='namespace':
				l=line_data[1].strip()
				r="%s\t%s\n" % (son_node,l)    
				dic[son_node] = l    
				f1.write(r)
    return dic
    
dic_go_namespace = id_namespaces()

#generate graph
def get_go_graph(filename) :
    G2=nx.DiGraph()
    a = []
    fr = open(filename,"r")
    #fw = open ('edges_error.txt','w')
    line = fr.readline()
    for line in fr :
        line_str = line.strip().split("\t")
        if len(line_str)==1 :
            G2.add_node(line_str[0])
        else :
            if dic_go_namespace[line_str[0]]==dic_go_namespace[line_str[1]]:
                G2.add_edge(line_str[0],line_str[1])
            else :
                #fw.write(line_str[0]+"\t"+line_str[1]+'\n')
                continue
           
    # delete nodes degree=0
    b = []
    for i in G2.nodes():
        if len(G2.predecessors(i))==0 and len(G2.successors(i))==0:
            a.append(i)
        if len(G2.predecessors(i))>0 and len(G2.successors(i))==0:
            b.append(i)
    print '移除孤立的点',len(a)
    #print '叶子节点',len(b)
    G2.remove_nodes_from(a)
    n=0
    sub_gragh = {}
    for i in nx.weakly_connected_component_subgraphs(G2):
        #print i.number_of_edges(),i.number_of_nodes()
        sub_gragh[n] = i
        n=n+1
   # print 'go_graph子网数',len(sub_gragh)
    #print [len(c) for c in sorted(nx.weakly_connected_components(G2), key=len, reverse=True)]
    return G2,b

def get_go_BP_graph(filename) :
    BP=nx.DiGraph()
    fr = open(filename,"r")
    for line in fr :
        line_str = line.strip().split("\t")
        node1, node2 = line_str
        BP.add_edge(node1,node2)
    fr.close()
    print 'BP:',BP.number_of_edges(),BP.number_of_nodes()
    return BP        
    
def get_gene_graph(filename) :
    G1=nx.Graph()
    fr = open(filename,"r")
    #fw = open ('edges_error.txt','w')
    for line in fr :
        line_str = line.strip().split("\t")
        G1.add_edge(line_str[0],line_str[1])
    n=0
    sub_gragh = {}
    for i in nx.connected_component_subgraphs(G1):
        sub_gragh[n] = i
        n=n+1
    #print 'gene_graph子网数',len(sub_gragh)
    #print [len(c) for c in sorted(nx.connected_component_subgraphs(G1), key=len, reverse=True)]

    return G1

#calculate depth
def get_term_depth(G2) :
    dic_term = {}
    set1 = set()
    set2 = set()
    set3 = set()
    length = nx.all_pairs_shortest_path_length(G2)
    err1=0
    for i in G2.nodes():
        try :
            dic_term[i] =length['GO:0003674'][i]
            set1.add(i)
        except :
            err1= err1+1  
    set_r=set()
    err2 = 0
    for i in G2.nodes():
        try :
            x=length['GO:0008150'][i]
            if not dic_term.has_key(i):
                dic_term[i] = x
                set2.add(i)
            else :
                set_r.add(i)
        except :
            err2= err2+1
    err3 = 0
    for i in G2.nodes():
        try :
            dic_term[i] =length['GO:0005575'][i]
            set3.add(i)
        except :
            err3= err3+1
    print 'dic_term is created:',len(set1),len(set2),len(set3)
    
    return dic_term

def get_BP_depth(BP) :
    
    dic_BPdepth = {}

    length = nx.all_pairs_shortest_path_length(BP)
    for node in BP.nodes() :
        x = length['GO:0008150'][node]
        #x = length['GO:0005575'][node]
        #x = length['GO:0003674'][node]
        dic_BPdepth[node] = x
    print 'dic_BPdepth is created:',len(dic_BPdepth)   
    return dic_BPdepth
    
#write depth to term_gene_depth.txt file    
def write_depth_to_file(dic_term):
    fr = open (r'G:\networkx\relationship\go_term_gene_one.txt','r')
    fw = open('term_gene_depth.txt','w')
    set_term_error = set()
    for line in fr :
        line_str = line.strip().split('\t')
        try:
            fw.write(str(dic_term[line_str[0].strip()])+'\t'+line_str[0].strip()+'\t'+line_str[1].strip()+'\n')
        except:
            fw.write('-1'+'\t'+line_str[0].strip()+'\t'+line_str[1].strip()+'\n')
            set_term_error.add(line_str[0])
    for i in set_term_error :
        print i
    print 'write success'
    fr.close()
    fw.close()
    
def draw_depth_Histogram(dic_term) :
    set4 =set()
    d = []
    for i in dic_term :
        set4.add(dic_term[i])
        d.append(dic_term[i])
    print '不同的层数:',len(set4)
    depth = []
    num = []
    for i in set4 :
        depth.append(i)
        num.append(d.count(i))
    #print i,d.count(i)
    plt.bar(depth,num,color='pink',width=1,align="center")
    plt.title('Depth_of_GoTree')
    plt.xlabel('depth')
    plt.ylabel('numbers')
#plt.legend('numbers',loc='upper')
    plt.xticks(depth)
    plt.grid(color='#95a5a6',linestyle='--', linewidth=1,axis='x',alpha=0.4)
    for i in range(0,len(num)):
        plt.text(depth[i]-0.25,num[i]+100,num[i])
    plt.show()

#the number of genes related to a term
def num_of_term_genes() :
    fr1 = open (r'G:\project1\term_gene\ccc.txt','r')
    fr2 = open (r'G:\project1\term_gene\fff.txt','r')
    fr3 = open (r'G:\project1\term_gene\ppp.txt','r')
    fw1 = open('ccc_num.txt','w')
    a = []
    dic= {}
    set1 = set()
    for line in fr1 :
        line_str = line.strip().split('\t')
        term = line_str[0].strip()
        #gene = line_str[1].strip()
        a.append(term)
        set1.add(term)
    for i in set1 :
        dic[i] = a.count(i)
    print len(dic)
    for i in dic :
        fw1.write(i+'\t'+str(dic[i])+'\n')
    fw1.close()
    fr1.close()
    fr2.close()
    fr3.close()
    
#not inherit term-gene
def gene_dic(G) :
    fr = open(r"G:\networkx\go_gene_dic.txt","r")
    #fw = open('term_geneNumber.txt','w')
    #fr1 = open (r'G:\networkx\go.txt','r')
    dic_gene = {}
    #dic_term = {}
    for line in fr:
        line_str = line.strip().split("*")
        term = line_str[0]
        gene = line_str[1]
        dic_gene[term] = []
        g = gene.strip().split('\t')
        for i in g :
            dic_gene[term].append(i.strip())
    for i in G.nodes() :
        if not dic_gene.has_key(i) :
            dic_gene[i] =[]
    return dic_gene
   
"""
G2 = get_go_graph(r"G:\networkx\relationship\go_go_one.txt")
dic_gene= gene_dic(G2)
n=0
sub_gragh = {}
for i in nx.weakly_connected_component_subgraphs(G2):
    sub_gragh[n] = i
    n=n+1
"""  
# inherit term-gene
tag = {}
def combine (node) :
    if  len(G2.successors(node))==0:
        return 
    else:
        for i in G2.successors(node):
            if tag.has_key(i) :
                dic_gene[node].extend(dic_gene[i])
            else:
                combine(i)
                dic_gene[node].extend(dic_gene[i])
                tag[i] =1

def gene_num() :
    fr = open(r'G:\project1\inherit_gene.txt','r')
    a = []
    dic = {}
    set1 = set()
    n = 0
    for line in fr :
        line_str = line.strip().split('\t')
        term = line_str[0].strip()
        gene_num = int(line_str[1].strip())
        dic[term] = gene_num
        set1.add(gene_num)
        a.append(gene_num)
        if gene_num == 50 :
            n=n+1
    print n
    print '不同的基因个数',len(set1)
    gene_num =[]
    term_num = []

    for i in set1:
        gene_num.append(i)
        term_num.append(a.count(i))
    return dic
        #print i,a.count(i),'个'
    #print sorted(term_num)
    #print sorted(gene_num)
    #plt.scatter(gene_num, term_num)


#找出注释基因为50的且没有继承关系的   
def find_term_by_geneNum(G,num1,num2,dic_term,filename)  :
    
    fr = open (filename,'r')
    b = []
    term = []
    dic_gene = {}
    line = fr.readline()
    for line in fr :
        line_str = line.strip().split('\t')
        termt = line_str[1].strip()
        geneNum= int(line_str[2].strip())
        if geneNum>=num1 and geneNum<=num2 :
            b.append(termt)
            dic_gene[termt] =[]
            for i in range(3,len(line_str)):
                dic_gene[termt].append(line_str[i])
    print '满足注释基因的term个数：',len(b)
    sub = G.subgraph(b)
    for i in nx.weakly_connected_component_subgraphs(sub):
        
        maxx = 0
        if i.number_of_nodes()>1:
            #print '构成子网的term',i.nodes()
            for t in i :
                depth = dic_term[t]
                if depth > maxx :
                    maxx = depth
                    node = t                  
            term.append(node)
            #print '在子网中选择靠近底部的term',node 
        else :
            term.extend(i)
    #print '差集(删掉的term)',list(set(b).difference(set(term))) 
    for i in list(set(b).difference(set(term))) :
        del dic_gene[i]
    relate_gene = []
    for i in dic_gene :
        relate_gene.extend(dic_gene[i])
    print '筛选过后的满足注释基因个数的term',len(term),'第一次涉及到的基因个数：',len(set(relate_gene))    
    return term,dic_gene
  
def find_subnet(dic_gene,G1,file_name) :
    fw = open(file_name,'w')
    subnets = {}
    set1 = set()
    for term in dic_gene :
        s_count = 0
        geneList = dic_gene[term]
        subnet = G1.subgraph(geneList)
        for s in nx.connected_component_subgraphs(subnet):     
            if s.number_of_nodes()>=5:    
                subnets[term+ '_' + str(s_count)] = s
                s_count = s_count + 1           
                for i in s.edges():
                    fw.write(term+'\t'+i[0]+'\t'+i[1]+'\n')
            if s.number_of_nodes()<5:
                for i in s.nodes() :
                    set1.add(i)
                    
       #print [len(c) for c in sorted(nx.connected_component_subgraphs(subnet), key=len, reverse=True)]
    print '筛选后的term涉及的基因子网数：',len(subnets),'最终涉及的基因个数',len(set1)
    fw.close()
   
    return subnets 
        
def  calculate_aggrement(filename1,filename2) :
    fr1 = open(filename1,'r')
    fr2 = open(filename2,'r')
    set1 = set()
    set2 = set()
    title1 = fr1.readline()
    title2 = fr2.readline()
    for line in fr1 :
        line_str = line.strip().split('\t')
        term = line_str[0].split('_')[0]
        set1.add(term)
    for line in fr2 :
        line_str = line.strip().split('\t')
        term = line_str[0].split('_')[0]
        set2.add(term)
    #for i in set1&set2:
        #print i
    print 'aggrement:',round(len(set1&set2)/float(len(set1|set2)),3)

def calculate_agreement_topN (filename1, filename2, N):
    fr1 = open(filename1,'r')
    fr2 = open(filename2,'r')
    line1 = fr1.readline()
    line2 = fr2.readline()
    set1= set()
    set2 = set()
    dic1 = {}
    dic2 = {}
    #fw = open("Conclusions_inter.csv","a")
    
    '''
    for line in fr1.readlines()[0:N] :
        line_str = line.strip().split('\t')
        set1.add(line_str[0].strip().split('_')[0])
    for line in fr2.readlines()[0:N] :
        line_str = line.strip().split('\t')
        set2.add(line_str[0].strip().split('_')[0])
    '''
    for line in fr1.readlines()[0:N] :
        line_str = line.strip().split(',')       
        term = line_str[0].strip().split('_')[0]
        if dic1.has_key(term):
            dic1[term] = dic1[term]+line_str[2:len(line_str)]
        else :
            dic1[term] = line_str[2:len(line_str)]
    for line in fr2.readlines()[0:N] :
        line_str = line.strip().split(',')
        term = line_str[0].strip().split('_')[0]
        if dic2.has_key(term):
            dic2[term] = dic2[term]+line_str[2:len(line_str)]
        else :
            dic2[term] = line_str[2:len(line_str)]
   
    intersect_subnet =  list(set(dic1.keys()).intersection(set(dic2.keys()))) #intersection of pathways names
    unino_subnet = list(set(dic1.keys()).union(set(dic2.keys()))) #union of pathways names
    g1 = []
    g2 = []
    dic_jiao = {}
    # caculate gene agreement 
    if len(intersect_subnet)==0 :
        agreement = 0
        #print 'top%d jiao_subnets:0' % (N)
        #fw.write('Top'+str(N)+','+'0'+'\n')
    else :
        for i in intersect_subnet :
            g1 = g1+dic1[i]     #file1 intersection genes
            g2 = g2+dic2[i]     #file2 intersection genes         
        intersect_gene  = set(g1) & set(g2)     #intersection genes between file1 and file2    
        union_gene = set(g1) | set(g2)   #union genes
        
        agreement = round(len(intersect_gene)/float(len(union_gene)),3)
        agreement1 = round(len(intersect_subnet)/float(len(unino_subnet)),3)          
        #print 'top%d jiao_subnets:%d jiao_Genes:%d gene_agree:%.3f' % (N,len(intersect_subnet),len(intersect_gene),agreement) 
        #fw.write('Top'+str(N)+','+str(len(intersect_subnet))+','+str(agreement1)+','+str(len(intersect_gene))+','+str(agreement)+'\n')
        sub = ','.join(intersect_subnet)
        #fw.write(','+sub+'\n')
        gene = ','.join(intersect_gene)
        #fw.write(','+gene+'\n')
    #print 'top',N,':',len(intersect_subnet),len(unino_subnet),len(intersect_gene),len(union_gene)     
    for i in intersect_subnet :
        dic_jiao [i] = list(set(dic1[i] + dic2[i]))
    fr1.close()
    fr2.close()
   
    '''
    agreement = len(set1 & set2)/float(len(set1 | set2))
    inters = set1&set2
    gene_set1 = set()
    gene_set2 = set()
    for line in fr1.readlines()[0:N]:
        line_str = line.strip().split('\t')
    #print sorted(set1&set2)
    #print sorted(set1|set2)
    print 'top%dIntersection: %d  subnet agreement:   %f' % (N,len(set1&set2),round(agreement,3))
    
    return set1&set2
    ''' 
    return dic_jiao

def chaji_gene(dic_geneN_geneI,dic1,dic2,term):
    
    a = dic1[term]
    b = dic2[term]

    chaji = list((set(a)|set(b)) - (set(a)&set(b)))
    print term,chaji
    for i in chaji :
        print i ,dic_geneN_geneI[i]
 
if __name__ =='__main__':

    dic_geneN_geneI = {}
    dic1 ={}
    dic2 = {}
    fr = open(r'G:\files\geneN_geneI.txt','r')
    for line in fr:
        gene,gid = line.strip().split('\t')
        dic_geneN_geneI[gid] = gene
    fr.close()
    dic1 =calculate_agreement_topN(r'G:\project1\BP\Pers_test.csv',r'G:\project1\BP\Haslett_test.csv',20)  
    dic2 = calculate_agreement_topN(r'G:\project1\BP\h_moregenes2\significant_subnetworks_test.csv',r'G:\project1\BP\p_moregenes2\significant_subnetworks_test.csv',20)
    print 'DMD'
    chaji_gene(dic_geneN_geneI,dic1,dic2,'GO:0019886')
    chaji_gene(dic_geneN_geneI,dic1,dic2,'GO:0051291')
    chaji_gene(dic_geneN_geneI,dic1,dic2,'GO:0060333')
    chaji_gene(dic_geneN_geneI,dic1,dic2,'GO:0060337')
    chaji_gene(dic_geneN_geneI,dic1,dic2,'GO:0032963')
    print 'Leuk'
    dic3 = calculate_agreement_topN(r'G:\project1\BP\o_arm_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\o_g_moregenes\significant_subnetworks_test.csv',20)
    dic4 = calculate_agreement_topN(r'G:\project1\BP\arm_moregenes2\significant_subnetworks_test.csv',r'G:\project1\BP\g_moregenes2\significant_subnetworks_test.csv',20)
    chaji_gene(dic_geneN_geneI,dic3,dic4,'GO:0048024')
    chaji_gene(dic_geneN_geneI,dic3,dic4,'GO:0031124')
    chaji_gene(dic_geneN_geneI,dic3,dic4,'GO:0046031')
    chaji_gene(dic_geneN_geneI,dic3,dic4,'GO:0006369')
    chaji_gene(dic_geneN_geneI,dic3,dic4,'GO:0046364')
    #chaji_gene(dic_geneN_geneI,dic1,dic2,'GO:0048024')
    chaji_gene(dic_geneN_geneI,dic3,dic4,'GO:0006090')

    '''
    #print 'original BP-control DMD!!!!!'  
    calculate_agreement_topN(r'G:\project1\BP\Pers_control.csv',r'G:\project1\BP\Haslett_control.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\Pers_control.csv',r'G:\project1\BP\Haslett_control.csv',20)    
    calculate_agreement_topN(r'G:\project1\BP\Pers_control.csv',r'G:\project1\BP\Haslett_control.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\Pers_control.csv',r'G:\project1\BP\Haslett_control.csv',40)

    #print 'original BP-test  DMD!!!!!!'
    calculate_agreement_topN(r'G:\project1\BP\Pers_test.csv',r'G:\project1\BP\Haslett_test.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\Pers_test.csv',r'G:\project1\BP\Haslett_test.csv',20)    
    calculate_agreement_topN(r'G:\project1\BP\Pers_test.csv',r'G:\project1\BP\Haslett_test.csv',30)         
    calculate_agreement_topN(r'G:\project1\BP\Pers_test.csv',r'G:\project1\BP\Haslett_test.csv',40)    

    
    '''
    '''
    print '1/log_density BP-control  log DMD !!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\Haslett_control_log_density.csv',r'G:\project1\BP\Pers_control_log_density.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\Haslett_control_log_density.csv',r'G:\project1\BP\Pers_control_log_density.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\Haslett_control_log_density.csv',r'G:\project1\BP\Pers_control_log_density.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\Haslett_control_log_density.csv',r'G:\project1\BP\Pers_control_log_density.csv',40)
    print '1/log_density BP-test log  DMD !!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\Haslett_test_log_density.csv',r'G:\project1\BP\Pers_test_log_density.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\Haslett_test_log_density.csv',r'G:\project1\BP\Pers_test_log_density.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\Haslett_test_log_density.csv',r'G:\project1\BP\Pers_test_log_density.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\Haslett_test_log_density.csv',r'G:\project1\BP\Pers_test_log_density.csv',40)       

    print '1/log_density BP-control more genes 1/2  DMD!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\p_moregenes\significant_subnetworks_control.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\p_moregenes\significant_subnetworks_control.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\p_moregenes\significant_subnetworks_control.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\p_moregenes\significant_subnetworks_control.csv',40)
    
    print '1/log_density BP-test more genes 1/2  DMD!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\p_moregenes\significant_subnetworks_test.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\p_moregenes\significant_subnetworks_test.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\p_moregenes\significant_subnetworks_test.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\p_moregenes\significant_subnetworks_test.csv',40)
    print '-----------------------------------------------------'
    '''
    '''
    #print '1/log_density BP-control more genes 2/3   DMD!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes2\significant_subnetworks_control.csv',r'G:\project1\BP\p_moregenes2\significant_subnetworks_control.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes2\significant_subnetworks_control.csv',r'G:\project1\BP\p_moregenes2\significant_subnetworks_control.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes2\significant_subnetworks_control.csv',r'G:\project1\BP\p_moregenes2\significant_subnetworks_control.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes2\significant_subnetworks_control.csv',r'G:\project1\BP\p_moregenes2\significant_subnetworks_control.csv',40)
    
    #print '1/log_density BP-test more genes 2/3   DMD!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes2\significant_subnetworks_test.csv',r'G:\project1\BP\p_moregenes2\significant_subnetworks_test.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes2\significant_subnetworks_test.csv',r'G:\project1\BP\p_moregenes2\significant_subnetworks_test.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes2\significant_subnetworks_test.csv',r'G:\project1\BP\p_moregenes2\significant_subnetworks_test.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes2\significant_subnetworks_test.csv',r'G:\project1\BP\p_moregenes2\significant_subnetworks_test.csv',40)
    '''


    #print '-----------------------------------------------------'
    '''
    print '1/log_density BP-control more genes 3/4 DMD!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes3\significant_subnetworks_control.csv',r'G:\project1\BP\p_moregenes3\significant_subnetworks_control.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes3\significant_subnetworks_control.csv',r'G:\project1\BP\p_moregenes3\significant_subnetworks_control.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes3\significant_subnetworks_control.csv',r'G:\project1\BP\p_moregenes3\significant_subnetworks_control.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes3\significant_subnetworks_control.csv',r'G:\project1\BP\p_moregenes3\significant_subnetworks_control.csv',40)
    
    print '1/log_density BP-test more genes3/4  DMD!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes3\significant_subnetworks_test.csv',r'G:\project1\BP\p_moregenes3\significant_subnetworks_test.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes3\significant_subnetworks_test.csv',r'G:\project1\BP\p_moregenes3\significant_subnetworks_test.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes3\significant_subnetworks_test.csv',r'G:\project1\BP\p_moregenes3\significant_subnetworks_test.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\h_moregenes3\significant_subnetworks_test.csv',r'G:\project1\BP\p_moregenes3\significant_subnetworks_test.csv',40)
    print '-----------------------------------------------------'
    '''
    '''
    #print '-----------------------------------------------------'
    #print 'original  control Leuk!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\o_arm_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\o_g_moregenes\significant_subnetworks_control.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\o_arm_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\o_g_moregenes\significant_subnetworks_control.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\o_arm_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\o_g_moregenes\significant_subnetworks_control.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\o_arm_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\o_g_moregenes\significant_subnetworks_control.csv',40)
    
    #print 'original  test Leuk !!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\o_arm_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\o_g_moregenes\significant_subnetworks_test.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\o_arm_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\o_g_moregenes\significant_subnetworks_test.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\o_arm_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\o_g_moregenes\significant_subnetworks_test.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\o_arm_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\o_g_moregenes\significant_subnetworks_test.csv',40)
    '''
    '''
    print '1/log_density BP-control log Leuk!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\log_arm_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\log_g_moregenes\significant_subnetworks_control.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\log_arm_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\log_g_moregenes\significant_subnetworks_control.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\log_arm_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\log_g_moregenes\significant_subnetworks_control.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\log_arm_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\log_g_moregenes\significant_subnetworks_control.csv',40)
    
    print '1/log_density BP-test log Leuk!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\log_arm_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\log_g_moregenes\significant_subnetworks_test.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\log_arm_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\log_g_moregenes\significant_subnetworks_test.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\log_arm_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\log_g_moregenes\significant_subnetworks_test.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\log_arm_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\log_g_moregenes\significant_subnetworks_test.csv',40)
    
    print '1/log_density BP-control more genes 1/2  Leuk!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\g_moregenes\significant_subnetworks_control.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\g_moregenes\significant_subnetworks_control.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\g_moregenes\significant_subnetworks_control.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\g_moregenes\significant_subnetworks_control.csv',40)
    
    print '1/log_density BP-test more genes 1/2  Leuk!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\g_moregenes\significant_subnetworks_test.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\g_moregenes\significant_subnetworks_test.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\g_moregenes\significant_subnetworks_test.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\g_moregenes\significant_subnetworks_test.csv',40)
    '''
    '''
   # print '1/log_density BP-control more genes 2/3  Leuk!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes2\significant_subnetworks_control.csv',r'G:\project1\BP\g_moregenes2\significant_subnetworks_control.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes2\significant_subnetworks_control.csv',r'G:\project1\BP\g_moregenes2\significant_subnetworks_control.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes2\significant_subnetworks_control.csv',r'G:\project1\BP\g_moregenes2\significant_subnetworks_control.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes2\significant_subnetworks_control.csv',r'G:\project1\BP\g_moregenes2\significant_subnetworks_control.csv',40)
    
    #print '1/log_density BP-test more genes 2/3  Leuk!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes2\significant_subnetworks_test.csv',r'G:\project1\BP\g_moregenes2\significant_subnetworks_test.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes2\significant_subnetworks_test.csv',r'G:\project1\BP\g_moregenes2\significant_subnetworks_test.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes2\significant_subnetworks_test.csv',r'G:\project1\BP\g_moregenes2\significant_subnetworks_test.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes2\significant_subnetworks_test.csv',r'G:\project1\BP\g_moregenes2\significant_subnetworks_test.csv',40)
    '''
    '''
    print '1/log_density BP-control more genes 3/4  Leuk!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes3\significant_subnetworks_control.csv',r'G:\project1\BP\g_moregenes3\significant_subnetworks_control.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes3\significant_subnetworks_control.csv',r'G:\project1\BP\g_moregenes3\significant_subnetworks_control.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes3\significant_subnetworks_control.csv',r'G:\project1\BP\g_moregenes3\significant_subnetworks_control.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes3\significant_subnetworks_control.csv',r'G:\project1\BP\g_moregenes3\significant_subnetworks_control.csv',40)
    
    print '1/log_density BP-test more genes 3/4   Leuk!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes3\significant_subnetworks_test.csv',r'G:\project1\BP\g_moregenes3\significant_subnetworks_test.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes3\significant_subnetworks_test.csv',r'G:\project1\BP\g_moregenes3\significant_subnetworks_test.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes3\significant_subnetworks_test.csv',r'G:\project1\BP\g_moregenes3\significant_subnetworks_test.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\arm_moregenes3\significant_subnetworks_test.csv',r'G:\project1\BP\g_moregenes3\significant_subnetworks_test.csv',40)

    print 'original---test---Pathway'
    print 'original pathway-control'
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\control_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\control_pers_density.csv',5)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\control_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\control_pers_density.csv',10)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\control_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\control_pers_density.csv',15)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\control_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\control_pers_density.csv',17)
    print 'original pathway-test'
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\test_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\test_pers_density.csv',5)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\test_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\test_pers_density.csv',10)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\test_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\test_pers_density.csv',15)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\test_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\test_pers_density.csv',17) 
    print 'log2 pathway-control'
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\control_haslett_density_log2.csv',r'G:\project1\pfsnet_results_top\density_log\control_pers_density_log2.csv',5)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\control_haslett_density_log2.csv',r'G:\project1\pfsnet_results_top\density_log\control_pers_density_log2.csv',10)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\control_haslett_density_log2.csv',r'G:\project1\pfsnet_results_top\density_log\control_pers_density_log2.csv',15)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\control_haslett_density_log2.csv',r'G:\project1\pfsnet_results_top\density_log\control_pers_density_log2.csv',17)
    print 'log2 pathway -test'
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\test_haslett_density_log2.csv',r'G:\project1\pfsnet_results_top\density_log\test_pers_density_log2.csv',5)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\test_haslett_density_log2.csv',r'G:\project1\pfsnet_results_top\density_log\test_pers_density_log2.csv',10)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\test_haslett_density_log2.csv',r'G:\project1\pfsnet_results_top\density_log\test_pers_density_log2.csv',15)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\test_haslett_density_log2.csv',r'G:\project1\pfsnet_results_top\density_log\test_pers_density_log2.csv',17)
    '''
    
    #calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\test_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\test_pers_density.csv',20)
    '''
    print 'dmd pathway more genes 2/3 control'    
    calculate_agreement_topN(r'G:\project1\path_has\significant_subnetworks_control.csv',r'G:\project1\path_pers\significant_subnetworks_control.csv',5)
    calculate_agreement_topN(r'G:\project1\path_has\significant_subnetworks_control.csv',r'G:\project1\path_pers\significant_subnetworks_control.csv',10)
    calculate_agreement_topN(r'G:\project1\path_has\significant_subnetworks_control.csv',r'G:\project1\path_pers\significant_subnetworks_control.csv',15)
    calculate_agreement_topN(r'G:\project1\path_has\significant_subnetworks_control.csv',r'G:\project1\path_pers\significant_subnetworks_control.csv',17)
    print 'dmd pathway more genes 2/3 test'
    calculate_agreement_topN(r'G:\project1\path_has\significant_subnetworks_test.csv',r'G:\project1\path_pers\significant_subnetworks_test.csv',5)
    calculate_agreement_topN(r'G:\project1\path_has\significant_subnetworks_test.csv',r'G:\project1\path_pers\significant_subnetworks_test.csv',10)
    calculate_agreement_topN(r'G:\project1\path_has\significant_subnetworks_test.csv',r'G:\project1\path_pers\significant_subnetworks_test.csv',15)
    calculate_agreement_topN(r'G:\project1\path_has\significant_subnetworks_test.csv',r'G:\project1\path_pers\significant_subnetworks_test.csv',17)
    
    print 'leuk pathway more genes 2/3 control'    
    calculate_agreement_topN(r'G:\project1\path_arm\significant_subnetworks_control.csv',r'G:\project1\path_gol\significant_subnetworks_control.csv',5)
    calculate_agreement_topN(r'G:\project1\path_arm\significant_subnetworks_control.csv',r'G:\project1\path_gol\significant_subnetworks_control.csv',10)
    calculate_agreement_topN(r'G:\project1\path_arm\significant_subnetworks_control.csv',r'G:\project1\path_gol\significant_subnetworks_control.csv',15)
    calculate_agreement_topN(r'G:\project1\path_arm\significant_subnetworks_control.csv',r'G:\project1\path_gol\significant_subnetworks_control.csv',17)
    print 'leuk pathway more genes 2/3 test'
    calculate_agreement_topN(r'G:\project1\path_arm\significant_subnetworks_test.csv',r'G:\project1\path_gol\significant_subnetworks_test.csv',5)
    calculate_agreement_topN(r'G:\project1\path_arm\significant_subnetworks_test.csv',r'G:\project1\path_gol\significant_subnetworks_test.csv',10)
    calculate_agreement_topN(r'G:\project1\path_arm\significant_subnetworks_test.csv',r'G:\project1\path_gol\significant_subnetworks_test.csv',15)
    calculate_agreement_topN(r'G:\project1\path_arm\significant_subnetworks_test.csv',r'G:\project1\path_gol\significant_subnetworks_test.csv',17)
   
    '''
    
    '''
    print 'control agreement level of Top N '    
    print 'original---'
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\control_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\control_pers_density.csv',5)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\control_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\control_pers_density.csv',10)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\control_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\control_pers_density.csv',15)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\control_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\control_pers_density.csv',17)
    #calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\control_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\control_pers_density.csv',20)
   
    print 'pow---'
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_pow\control_haslett_density_pow2.csv',r'G:\project1\pfsnet_results_top\density_pow\control_pers_density_pow2.csv',5)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_pow\control_haslett_density_pow2.csv',r'G:\project1\pfsnet_results_top\density_pow\control_pers_density_pow2.csv',10)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_pow\control_haslett_density_pow2.csv',r'G:\project1\pfsnet_results_top\density_pow\control_pers_density_pow2.csv',15)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_pow\control_haslett_density_pow2.csv',r'G:\project1\pfsnet_results_top\density_pow\control_pers_density_pow2.csv',17)
    #calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_pow\control_haslett_density_pow2.csv',r'G:\project1\pfsnet_results_top\density_pow\control_pers_density_pow2.csv',20)
    print 'exp---'
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_exp\control_haslett_density_exp.csv',r'G:\project1\pfsnet_results_top\density_exp\control_pers_density_exp.csv',5)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_exp\control_haslett_density_exp.csv',r'G:\project1\pfsnet_results_top\density_exp\control_pers_density_exp.csv',10)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_exp\control_haslett_density_exp.csv',r'G:\project1\pfsnet_results_top\density_exp\control_pers_density_exp.csv',15)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_exp\control_haslett_density_exp.csv',r'G:\project1\pfsnet_results_top\density_exp\control_pers_density_exp.csv',17)
    #calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_exp\control_haslett_density_exp.csv',r'G:\project1\pfsnet_results_top\density_exp\control_pers_density_exp.csv',20)
    print 'log---'
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\control_haslett_density_log.csv',r'G:\project1\pfsnet_results_top\density_log\control_pers_density_log.csv',5)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\control_haslett_density_log.csv',r'G:\project1\pfsnet_results_top\density_log\control_pers_density_log.csv',10)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\control_haslett_density_log.csv',r'G:\project1\pfsnet_results_top\density_log\control_pers_density_log.csv',15)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\control_haslett_density_log.csv',r'G:\project1\pfsnet_results_top\density_log\control_pers_density_log.csv',17)
    #calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\control_haslett_density_log.csv',r'G:\project1\pfsnet_results_top\density_log\control_pers_density_log.csv',20)
    print 'density---'
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density\control_haslett_density.csv',r'G:\project1\pfsnet_results_top\density\control_pers_density.csv',5)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density\control_haslett_density.csv',r'G:\project1\pfsnet_results_top\density\control_pers_density.csv',10)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density\control_haslett_density.csv',r'G:\project1\pfsnet_results_top\density\control_pers_density.csv',15)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density\control_haslett_density.csv',r'G:\project1\pfsnet_results_top\density\control_pers_density.csv',17)
    #calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density\control_haslett_density.csv',r'G:\project1\pfsnet_results_top\density\control_pers_density.csv',20)
    '''

    '''
    print 'GO original control ALL!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\o_a_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\o_m_moregenes\significant_subnetworks_control.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\o_a_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\o_m_moregenes\significant_subnetworks_control.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\o_a_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\o_m_moregenes\significant_subnetworks_control.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\o_a_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\o_m_moregenes\significant_subnetworks_control.csv',40)
    
    print 'original  test ALL !!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\o_a_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\o_m_moregenes\significant_subnetworks_test.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\o_a_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\o_m_moregenes\significant_subnetworks_test.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\o_a_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\o_m_moregenes\significant_subnetworks_test.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\o_a_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\o_m_moregenes\significant_subnetworks_test.csv',40)
    
    
    print '1/log_density BP-control log  ALL!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\log_a_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\log_m_moregenes\significant_subnetworks_control.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\log_a_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\log_m_moregenes\significant_subnetworks_control.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\log_a_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\log_m_moregenes\significant_subnetworks_control.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\log_a_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\log_m_moregenes\significant_subnetworks_control.csv',40)
    
    print '1/log_density BP-test log  ALL !!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\log_a_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\log_m_moregenes\significant_subnetworks_test.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\log_a_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\log_m_moregenes\significant_subnetworks_test.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\log_a_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\log_m_moregenes\significant_subnetworks_test.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\log_a_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\log_m_moregenes\significant_subnetworks_test.csv',40)

    print '1/log_density BP-control more genes  ALL!!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\a_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\m_moregenes\significant_subnetworks_control.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\a_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\m_moregenes\significant_subnetworks_control.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\a_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\m_moregenes\significant_subnetworks_control.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\a_moregenes\significant_subnetworks_control.csv',r'G:\project1\BP\m_moregenes\significant_subnetworks_control.csv',40)
    
    print '1/log_density BP-test more genes   ALL !!!!!!!!!'    
    calculate_agreement_topN(r'G:\project1\BP\a_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\m_moregenes\significant_subnetworks_test.csv',10)
    calculate_agreement_topN(r'G:\project1\BP\a_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\m_moregenes\significant_subnetworks_test.csv',20)
    calculate_agreement_topN(r'G:\project1\BP\a_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\m_moregenes\significant_subnetworks_test.csv',30)
    calculate_agreement_topN(r'G:\project1\BP\a_moregenes\significant_subnetworks_test.csv',r'G:\project1\BP\m_moregenes\significant_subnetworks_test.csv',40)
    '''
    '''
    print 'test agreement level of Top N '
    print 'original---'
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\test_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\test_pers_density.csv',5)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\test_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\test_pers_density.csv',10)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\test_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\test_pers_density.csv',15)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\test_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\test_pers_density.csv',17)
    #calculate_agreement_topN(r'G:\project1\pfsnet_results_top\original\test_haslett_density.csv',r'G:\project1\pfsnet_results_top\original\test_pers_density.csv',20)
    print 'pow---'
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_pow\test_haslett_density_pow2.csv',r'G:\project1\pfsnet_results_top\density_pow\test_pers_density_pow2.csv',5)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_pow\test_haslett_density_pow2.csv',r'G:\project1\pfsnet_results_top\density_pow\test_pers_density_pow2.csv',10)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_pow\test_haslett_density_pow2.csv',r'G:\project1\pfsnet_results_top\density_pow\test_pers_density_pow2.csv',15)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_pow\test_haslett_density_pow2.csv',r'G:\project1\pfsnet_results_top\density_pow\test_pers_density_pow2.csv',17)
    #calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_pow\test_haslett_density_pow2.csv',r'G:\project1\pfsnet_results_top\density_pow\test_pers_density_pow2.csv',20)
    print 'exp---'
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_exp\test_haslett_density_exp.csv',r'G:\project1\pfsnet_results_top\density_exp\test_pers_density_exp.csv',5)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_exp\test_haslett_density_exp.csv',r'G:\project1\pfsnet_results_top\density_exp\test_pers_density_exp.csv',10)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_exp\test_haslett_density_exp.csv',r'G:\project1\pfsnet_results_top\density_exp\test_pers_density_exp.csv',15)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_exp\test_haslett_density_exp.csv',r'G:\project1\pfsnet_results_top\density_exp\test_pers_density_exp.csv',17)
    #calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_exp\test_haslett_density_exp.csv',r'G:\project1\pfsnet_results_top\density_exp\test_pers_density_exp.csv',20)
    print 'log---'
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\test_haslett_density_log.csv',r'G:\project1\pfsnet_results_top\density_log\test_pers_density_log.csv',5)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\test_haslett_density_log.csv',r'G:\project1\pfsnet_results_top\density_log\test_pers_density_log.csv',10)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\test_haslett_density_log.csv',r'G:\project1\pfsnet_results_top\density_log\test_pers_density_log.csv',15)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\test_haslett_density_log.csv',r'G:\project1\pfsnet_results_top\density_log\test_pers_density_log.csv',17)
    #calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density_log\test_haslett_density_log.csv',r'G:\project1\pfsnet_results_top\density_log\test_pers_density_log.csv',20)
    print 'density---'
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density\test_haslett_density.csv',r'G:\project1\pfsnet_results_top\density\test_pers_density.csv',5)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density\test_haslett_density.csv',r'G:\project1\pfsnet_results_top\density\test_pers_density.csv',10)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density\test_haslett_density.csv',r'G:\project1\pfsnet_results_top\density\test_pers_density.csv',15)
    calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density\test_haslett_density.csv',r'G:\project1\pfsnet_results_top\density\test_pers_density.csv',17)
    #calculate_agreement_topN(r'G:\project1\pfsnet_results_top\density\test_haslett_density.csv',r'G:\project1\pfsnet_results_top\density\test_pers_density.csv',20)    
    '''
    
    #BP = get_go_BP_graph(r'G:\project1\GO_tree\BP_edges.txt')
    #dic_BPdepth = get_BP_depth(BP)       
    #term ,dic_gene= find_term_by_geneNum(BP,50,100,dic_BPdepth,r'BP_Depth_NumInheritGene_Id.txt')
    #G1 = get_gene_graph(r'G:\project1\GoFunction_GeneId.txt')
    #subnets = find_subnet(dic_gene,G1,'BP50_100.txt')
    
    '''
    #+density
    print 'density agreement:(control/test/all)'    
    calculate_aggrement(r'G:\project1\density_p_value\pfsnet_results_haslet\significant_subnetworks_control.csv',r'G:\project1\density_p_value\pfsnet_results_pers\significant_subnetworks_control.csv')
    calculate_aggrement(r'G:\project1\density_p_value\pfsnet_results_pers\significant_subnetworks_test.csv',r'G:\project1\density_p_value\pfsnet_results_haslet\significant_subnetworks_test.csv')    
    calculate_aggrement(r'G:\project1\density_p_value\pfsnet_results_pers\sub.csv',r'G:\project1\density_p_value\pfsnet_results_haslet\sub.csv')   
    #log density
    print 'density log agreement:(control/test/all)'    
    calculate_aggrement(r'G:\project1\density_p_value_log\pfsnet_results_haslet\significant_subnetworks_control.csv',r'G:\project1\density_p_value_log\pfsnet_results_pers\significant_subnetworks_control.csv')
    calculate_aggrement(r'G:\project1\density_p_value_log\pfsnet_results_pers\significant_subnetworks_test.csv',r'G:\project1\density_p_value_log\pfsnet_results_haslet\significant_subnetworks_test.csv')    
    calculate_aggrement(r'G:\project1\density_p_value_log\pfsnet_results_pers\sub.csv',r'G:\project1\density_p_value_log\pfsnet_results_haslet\sub.csv')   
    print 'density log1 agreement:(control/test/all)'    
    calculate_aggrement(r'G:\project1\density_p_value_log1\pfsnet_results_haslet\significant_subnetworks_control.csv',r'G:\project1\density_p_value_log1\pfsnet_results_pers\significant_subnetworks_control.csv')
    calculate_aggrement(r'G:\project1\density_p_value_log1\pfsnet_results_pers\significant_subnetworks_test.csv',r'G:\project1\density_p_value_log1\pfsnet_results_haslet\significant_subnetworks_test.csv')    
    calculate_aggrement(r'G:\project1\density_p_value_log1\pfsnet_results_pers\sub.csv',r'G:\project1\density_p_value_log1\pfsnet_results_haslet\sub.csv')   
    print 'density log2 agreement:(control/test/all)'    
    calculate_aggrement(r'G:\project1\density_p_value_log2\pfsnet_results_haslet\significant_subnetworks_control.csv',r'G:\project1\density_p_value_log2\pfsnet_results_pers\significant_subnetworks_control.csv')
    calculate_aggrement(r'G:\project1\density_p_value_log2\pfsnet_results_pers\significant_subnetworks_test.csv',r'G:\project1\density_p_value_log2\pfsnet_results_haslet\significant_subnetworks_test.csv')    
    calculate_aggrement(r'G:\project1\density_p_value_log2\pfsnet_results_pers\sub.csv',r'G:\project1\density_p_value_log2\pfsnet_results_haslet\sub.csv')   
    print 'density log3 agreement:(control/test/all)'    
    calculate_aggrement(r'G:\project1\density_p_value_log3\pfsnet_results_haslet\significant_subnetworks_control.csv',r'G:\project1\density_p_value_log3\pfsnet_results_pers\significant_subnetworks_control.csv')
    calculate_aggrement(r'G:\project1\density_p_value_log3\pfsnet_results_pers\significant_subnetworks_test.csv',r'G:\project1\density_p_value_log3\pfsnet_results_haslet\significant_subnetworks_test.csv')    
    calculate_aggrement(r'G:\project1\density_p_value_log3\pfsnet_results_pers\sub.csv',r'G:\project1\density_p_value_log3\pfsnet_results_haslet\sub.csv')   
    
    print 'density exp agreement:(control/test/all)'    
    calculate_aggrement(r'G:\project1\density_p_value_exp\pfsnet_results_haslet\significant_subnetworks_control.csv',r'G:\project1\density_p_value_exp\pfsnet_results_pers\significant_subnetworks_control.csv')
    calculate_aggrement(r'G:\project1\density_p_value_exp\pfsnet_results_pers\significant_subnetworks_test.csv',r'G:\project1\density_p_value_exp\pfsnet_results_haslet\significant_subnetworks_test.csv')    
    calculate_aggrement(r'G:\project1\density_p_value_exp\pfsnet_results_pers\sub.csv',r'G:\project1\density_p_value_exp\pfsnet_results_haslet\sub.csv')   
   
    
    print 'density pow agreement:(control/test/all)'    
    calculate_aggrement(r'G:\project1\density_p_value_pow\pfsnet_results_haslet\significant_subnetworks_control.csv',r'G:\project1\density_p_value_pow\pfsnet_results_pers\significant_subnetworks_control.csv')
    calculate_aggrement(r'G:\project1\density_p_value_pow\pfsnet_results_pers\significant_subnetworks_test.csv',r'G:\project1\density_p_value_pow\pfsnet_results_haslet\significant_subnetworks_test.csv')    
    calculate_aggrement(r'G:\project1\density_p_value_pow\pfsnet_results_pers\sub.csv',r'G:\project1\density_p_value_pow\pfsnet_results_haslet\sub.csv')    

    #print 'original agreement:(control/test/all)'
    #calculate_aggrement(r'G:\project1\pfsnet_result_normal\pfsnet_results_haslet\significant_subnetworks_control.csv',r'G:\project1\pfsnet_result_normal\pfsnet_results_pers\significant_subnetworks_control.csv')
    #calculate_aggrement(r'G:\project1\pfsnet_result_normal\pfsnet_results_haslet\significant_subnetworks_test.csv',r'G:\project1\pfsnet_result_normal\pfsnet_results_pers\significant_subnetworks_test.csv')    
    #calculate_aggrement(r'G:\project1\pfsnet_result_normal\pfsnet_results_haslet\sub.csv',r'G:\project1\pfsnet_result_normal\pfsnet_results_pers\sub.csv') 

    print 'original1 agreement:(control/test/all)'
    calculate_aggrement(r'G:\project1\pfsnet_result_normal1\pfsnet_results_haslet\significant_subnetworks_control.csv',r'G:\project1\pfsnet_result_normal1\pfsnet_results_pers\significant_subnetworks_control.csv')
    calculate_aggrement(r'G:\project1\pfsnet_result_normal1\pfsnet_results_haslet\significant_subnetworks_test.csv',r'G:\project1\pfsnet_result_normal1\pfsnet_results_pers\significant_subnetworks_test.csv')    
    calculate_aggrement(r'G:\project1\pfsnet_result_normal1\pfsnet_results_haslet\sub.csv',r'G:\project1\pfsnet_result_normal1\pfsnet_results_pers\sub.csv') 

    print 'original2 agreement:(control/test/all)'
    calculate_aggrement(r'G:\project1\pfsnet_result_normal2\pfsnet_results_haslet\significant_subnetworks_control.csv',r'G:\project1\pfsnet_result_normal2\pfsnet_results_pers\significant_subnetworks_control.csv')
    calculate_aggrement(r'G:\project1\pfsnet_result_normal2\pfsnet_results_haslet\significant_subnetworks_test.csv',r'G:\project1\pfsnet_result_normal2\pfsnet_results_pers\significant_subnetworks_test.csv')    
    calculate_aggrement(r'G:\project1\pfsnet_result_normal2\pfsnet_results_haslet\sub.csv',r'G:\project1\pfsnet_result_normal2\pfsnet_results_pers\sub.csv') 
    print 'original3 agreement:(control/test/all)'
    calculate_aggrement(r'G:\project1\pfsnet_result_normal3\pfsnet_results_haslet\significant_subnetworks_control.csv',r'G:\project1\pfsnet_result_normal3\pfsnet_results_pers\significant_subnetworks_control.csv')
    calculate_aggrement(r'G:\project1\pfsnet_result_normal3\pfsnet_results_haslet\significant_subnetworks_test.csv',r'G:\project1\pfsnet_result_normal3\pfsnet_results_pers\significant_subnetworks_test.csv')    
    calculate_aggrement(r'G:\project1\pfsnet_result_normal3\pfsnet_results_haslet\sub.csv',r'G:\project1\pfsnet_result_normal3\pfsnet_results_pers\sub.csv') 


    #print 'BP_pow agreement:(control/test/all)'
    #calculate_aggrement(r'G:\project1\BP\pfsnet_results_BP_pow_haslet\significant_subnetworks_control.csv',r'G:\project1\BP\pfsnet_results_BP_pow_pers\significant_subnetworks_control.csv')
    #calculate_aggrement(r'G:\project1\BP\pfsnet_results_BP_pow_haslet\significant_subnetworks_test.csv',r'G:\project1\BP\pfsnet_results_BP_pow_pers\significant_subnetworks_test.csv')
    '''
    '''
   
    G,b= get_go_graph(r"G:\networkx\relationship\go_go_one.txt") 
    dic_term = get_term_depth(G)  
    term ,dic_gene= find_term_by_geneNum(G,50,100,dic_term,r'Term_Depth_NumInheritGene_Id.txt')
    G1 = get_gene_graph(r'G:\project1\GoFunction_GeneId.txt')
    subnets = find_subnet(dic_gene,G1,'g50_100.txt')
    '''
    


    
    
"""
#注释基因为0的term分布情况统计
dic = gene_num()
n=0
m=0
k=0
c = set()
G2,b= get_go_graph(r"G:\networkx\relationship\go_go_one.txt")
for x in b:     
    for y in G2.predecessors(x):
        c.add(y)     
print len(b) 
for i in b :
    if dic[i]==0 :
        n=n+1
for i in c :
    if dic[i]==0 :
        m=m+1

print '叶子节点中注释基因=00的term个数',n
print '叶子节点的父亲结点中 注释基因=0的term个数',m
"""        
                        

"""

combine('GO:0008150') 
fw = open('inherit_gene3.txt','w') 
for i in sub_gragh[0].nodes():
    fw.write(i+'\t'+str(len(set(dic_gene[i])))+'\t') 
    for gene in set(dic_gene[i]):
        fw.write(gene+'\t')
    fw.write('\n')
fw.close()
"""

    
"""  
#Execute the statement  
G2 = get_go_graph(r"G:\networkx\relationship\go_go_one.txt")   
dic_term = get_term_depth(G2)
write_depth_to_file(dic_term)
draw_depth_Histogram(dic_term)
print G2.number_of_nodes(),'结点数'
print G2.number_of_edges(),'边数'
gene_inherit(G2)  

"""


"""
fr = open(r'G:\project1\inherit_gene.txt','r')
fw = open('Term_Depth_NumInheritGene.txt','w')  
G2 ,b= get_go_graph(r"G:\networkx\relationship\go_go_one.txt")  
dic_term = get_term_depth(G2) 
n=0
for line in fr :
    line_str = line.strip().split('\t')
    inherit_num_gene = int(line_str[1].strip())
    if inherit_num_gene >0: 
        fw.write(str(dic_term[line_str[0].strip()])+'\t'+line)
    else :
        n=n+1
fw.close()
print '注释基因为0的term数',n
"""
"""
输出两点之间的路径
paths = list(nx.shortest_simple_paths(G, 'GO:0005575', 'GO:0045271'))
for i in paths :
    print len(i)
    if len(i)==3 or len(i)==4 :
        print i
"""

