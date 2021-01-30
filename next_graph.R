library(igraph)

#1)导入边数据和节点数据
edges <- read.table('C:\\Users\\CHARLES\\OneDrive\\分析流程\\circrna\\DAY3_vs_WT5\\microrna_target\\mmu.txt', header=T, sep='\t')
#导入边数据，里面可以包含每个边的频次数据或权重
vertices <- read.table('C:\\Users\\CHARLES\\OneDrive\\分析流程\\circrna\\DAY3_vs_WT5\\microrna_target\\vertices.txt', header=T, sep='\t') #导入节点数据，可以包含属性数据，如分类

#2)导入数据后，要转化成图数据才能用R作图，不同数据格式用不同方式
#graph_from_literal 通过文字创建，graph_from_edgelist通过边列表（矩阵）创建，graph_from_adj_list通过邻接列表（矩阵）创建
#graph_from_adjacency_matrix 通过邻接矩阵（所有点的纵横矩阵）创建，graph_from_incidence_matrix通过关联矩阵（两组内部独立点的矩阵）创建
#graph_from_data_frame通过数据框创建，详细介绍可参见《igraph manual》http://igraph.org/r/doc/igraph.pdf
#这里数据格式是数据框，所以用graph_from_data_frame
graph <- graph_from_data_frame(edges, directed = TRUE, vertices=vertices) 
#directed = TRUE表示有方向,如果不需要点数据，可以设置vertices=NULL

#3)转换完成后，有两种生成方式，一是直接plot,参数放里面；二是通过修改图的方式设置参数，然后plot

#生成方式1：
plot(graph,  
	layout=layout.fruchterman.reingold,  #layout.fruchterman.reingold表示弹簧式发散的布局，
	#其他还有环形布局layout.circle，分层布局layout.reingold.tilford，中心向外发散layout.reingold.tilford(graph,circular=T) ，核心布局layout_as_star，大型网络可视化layout_with_drl
	vertex.size=5,						#节点大小  
	vertex.shape='circle',		#节点不带边框none,,圆形边框circle,方块形rectangle  
	vertex.color="yellow",		#设置颜色，其他如red,blue,cyan,yellow等
	vertex.label=NULL,				#NULL表示不设置，为默认状态，即默认显示数据中点的名称，可以是中文。如果是NA则表示不显示任何点信息	 
	vertex.label.cex=0.8,			#节点字体大小
	vertex.label.color='black',	#节点字体颜色,red  
	vertex.label.dist=0.4,		#标签和节点位置错开
	edge.arrow.size=0.3,			#连线的箭头的大小,若为0即为无向图，当然有些数据格式不支持有向图  
	edge.width = 0.5,					#连接线宽度
	edge.label=NA,						#不显示连接线标签，默认为频次
	edge.color="black")				#连线颜色

#生成方式2：
library(igraph)
edges <- read.table('C:\\Users\\CHARLES\\OneDrive\\分析流程\\circrna\\DAY7_vs_WT5\\microrna_target\\mmu.txt', header=T, sep='\t')
vertices <- read.table('C:\\Users\\CHARLES\\OneDrive\\分析流程\\circrna\\DAY7_vs_WT5\\microrna_target\\vertices.txt', header=T, sep='\t') #导入节点数据，可以包含属性数据，如分类
graph <- graph_from_data_frame(edges, directed = TRUE, vertices=vertices)
set.seed(500) #生成随机数，这样图的布局就会可重复，而不是每次生成的时候都变
l <- layout.fruchterman.reingold(graph) #设置图的布局方式为弹簧式发散的布局

#具体修改过程
#V(graph)$size <- degree(graph) * 2 + 5  #节点大小与点中心度成正比，中心度即与该点相连的点的总数
V(graph)$size <- 5;
#colrs <- c("red", "green")

#V(graph)$color <- V(graph)$color #根据类型设置颜色,按照类型分组
V(graph)$label <- NA
V(graph)$label.color <- "black" #设置节点标记的颜色
E(graph)$color="black"
#E(graph)$width <- E(graph)$fre #根据频次列设置边宽度
#E(graph)$label <- E(graph)$fre #根据频次列设置边标签
E(graph)$arrow.size=0.3 #设置箭头大小
#生成图
plot(graph, layout=l)
