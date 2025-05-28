#' Use the layout in the sne package to calculate the visual layout of the network
#'
#' @param cor Correlation matrix
#' @param method method to culculate Degree of modularity
#' @examples
#' data(igraph)
#' res = ZiPiPlot(igraph = igraph,method = "cluster_fast_greedy")
#' res[[1]]
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export



ZiPiPlot = function(igraph = igraph,method = "cluster_fast_greedy"){

  if (method == "cluster_walktrap" ) {
    fc <- igraph::cluster_walktrap(igraph,weights =  abs(igraph::E(igraph)$weight))
  }
  if (method == "cluster_edge_betweenness" ) {
    fc <- igraph::cluster_edge_betweenness(igraph,weights =  abs(igraph::E(igraph)$weight))
  }
  if (method == "cluster_fast_greedy" ) {
    fc <- igraph::cluster_fast_greedy(igraph,weights =  abs(igraph::E(igraph)$weight))
  }
  if (method == "cluster_spinglass" ) {
    fc <- igraph::cluster_spinglass(igraph,weights =  abs(igraph::E(igraph)$weight))
  }

  modularity <- igraph::modularity(igraph,igraph::membership(fc))
  # 模块化程度
  # 按照模块为节点配色，这里我们可以加入到nodes中
  comps <- igraph::membership(fc)

  # comps
  igraph::V(igraph)$module <- as.character(comps)

  taxa.roles <- module.roles(igraph)

  taxa.roles$label = row.names(taxa.roles)
  for (i in 1:nrow(taxa.roles))if(taxa.roles[i,3]> 0.62|taxa.roles[i,1]> 2.5) {
    taxa.roles[i,5]=taxa.roles[i,5]
  }else{
    taxa.roles[i,5]= ""
  }
  taxa.roles$role_7 = taxa.roles$roles

  taxa.roles <- na.omit(taxa.roles)   # remove NA values
  taxa.roles[which(taxa.roles$z < 2.5 & taxa.roles$p < 0.62),'roles'] <- 'Peripherals'
  taxa.roles[which(taxa.roles$z < 2.5 & taxa.roles$p >= 0.62),'roles'] <- 'Connectors'
  taxa.roles[which(taxa.roles$z >= 2.5 & taxa.roles$p < 0.62),'roles'] <- 'Module hubs'
  taxa.roles[which(taxa.roles$z >= 2.5 & taxa.roles$p >= 0.62),'roles'] <- 'Network hubs'



  p <- plot_roles2(taxa.roles) +
    ggrepel::geom_text_repel(data = taxa.roles,
                             aes(x = p, y = z, color = module,label=taxa.roles$label),size=4)#
  p
  #geom_text(data = taxa.roles, aes(x = p, y = z, color = module,label=taxa.roles$label),size=4)
  # print(p)

  return(list(# 定义优化的ZIPI函数----
ZiPi_v1 <- function(netw, modules=F) {
  if (!modules){
    community <- cluster_louvain(netw)
    modules <- membership(community)
  }
  # 提取网络节点名称
  names <- V(netw)$name
  
  # 提取邻接矩阵
  # adj_matrix <- as_adjacency_matrix(netw, sparse = FALSE)
  adj_matrix <- as_adj(netw, sparse = FALSE)
  
  # 总连接数和模块连接数计算
  total_connections <- rowSums(adj_matrix > 0)
  module_connections <- rowSums(adj_matrix * outer(modules, modules, "=="))
  
  # 计算 P 值
  unique_modules <- unique(modules)
  KitKi <- matrix(0, nrow = length(names), ncol = length(unique_modules))
  for (j in seq_along(unique_modules)) {
    module_mask <- (modules == unique_modules[j])
    KitKi[, j] <- rowSums(adj_matrix[, module_mask]) / total_connections
  }
  P <- 1 - rowSums(KitKi^2)
  
  # 计算模块的均值和标准差
  module_stats <- aggregate(module_connections, by = list(module = modules), FUN = function(x) c(mean = mean(x), sd = sd(x)))
  meanZ <- module_stats$x[, "mean"]
  sdZ <- module_stats$x[, "sd"]
  
  # 计算 Z 值
  Z <- (module_connections - meanZ[modules]) / sdZ[modules]
  
  # 返回结果
  return(data.frame(names, module = modules, module_connections, total_connections, z=Z, p=P))
}

ZiPi_v2 <- function(netw, modules = NULL) {
  # 如果未提供模块划分，使用 Louvain 方法进行模块检测
  if (is.null(modules)) {
    # 复制网络避免修改原数据
    netw_pos <- netw

    # 由于模块聚类方法不允许负值存在，采用线性平移转化为正值（保留方向性，需确保最小权重=0）
    min_weight <- min(E(netw)$weight)
    E(netw_pos)$weight <- E(netw)$weight - min_weight

    community <- cluster_louvain(netw_pos)
    modules <- membership(community)
  }
  

  # 提取网络节点名称和邻接矩阵
  node_names <- V(netw)$name
  adj_matrix <- as_adj(netw, sparse = FALSE)
  
  # 计算每个节点的总连接数
  total_connections <- rowSums(adj_matrix > 0)
  
  # 剔除模块大小为 1 的模块
  module_sizes <- table(modules)
  valid_modules <- names(module_sizes[module_sizes > 1])
  valid_mask <- modules %in% valid_modules
  
  # 筛选有效的节点和模块
  modules <- modules[valid_mask]
  adj_matrix <- adj_matrix[valid_mask, valid_mask, drop = FALSE]
  total_connections <- total_connections[valid_mask]
  node_names <- node_names[valid_mask]
  
  # 重新计算模块连接数
  module_connections <- rowSums(adj_matrix * outer(modules, modules, "=="))
  
  # 计算 P 值（模块连接分布均匀性）
  unique_modules <- unique(modules)
  KitKi <- matrix(0, nrow = length(node_names), ncol = length(unique_modules))
  for (j in seq_along(unique_modules)) {
    module_mask <- (modules == unique_modules[j])
    KitKi[, j] <- rowSums(adj_matrix[, module_mask]) / total_connections
  }
  P <- 1 - rowSums(KitKi^2)
  
  # 计算模块的均值和标准差
  module_stats <- aggregate(module_connections, by = list(module = modules), 
                            FUN = function(x) c(mean = mean(x), sd = sd(x)))
  meanZ <- module_stats$x[, "mean"]
  sdZ <- module_stats$x[, "sd"]
  
  # 计算 Z 值（节点在模块中的连通性）
  Z <- (module_connections - meanZ[modules]) / sdZ[modules]
  
  # 返回结果
  return(data.frame(
    name = node_names, 
    module = modules, 
    module_connections = module_connections, 
    total_connections = total_connections, 
    z = Z, 
    p = P
  ))
}


# 直接从zipi文件开始进行分析----
# zi_pi_BF3 <- read.csv("./psBF_AVS0913/BFZiPi.csv")
# zi_pi_WF3 <- read.csv("./psWF_AVS0913/WFZiPi.csv")
# zi_pi_BF2 <- read.csv("./psBFnoEnv_AVS0913/BFZiPi.csv")
# zi_pi_WF2 <- read.csv("./psWFnoEnv_AVS0913/WFZiPi.csv")
# zi_pi_merge3 <- read.csv("./psMerge_AVS0913/oneZiPi.csv")
# zi_pi_merge2 <- read.csv("./psMergenoEnv_AVS0913/oneZiPi.csv")
# zi_pi <- zi_pi_merge3
# head(zi_pi_BF2)
# str(igraph)
# zi_pi <- ZiPi_v2(netw=igraph)

# #可自定义module顺序，这里随意举例
# Lv <- as.character(unique(zi_pi$module))
# zi_pi$module <- factor(zi_pi$module, levels = Lv)
# colnames(zi_pi)[colnames(zi_pi) %in% c("Z","P")] <- c("z","p")
# zipi_redraw(zi_pi)

# 定义zipi_redraw绘图函数----
zipi_redraw <- function(zi_pi = zipi, levels = Lv, lable = F) #注意此处根据roles上色，即zipi表格中一定要有roles一列
{
  zi_pi$module <- factor(zi_pi$module, levels = Lv)
  ggplot(zi_pi, aes(p, z)) +
    # 注意，这里的roles如果为数字，aes上色会将其变为连续型变量，应该指定其为离散型
    geom_point(aes(color = module), alpha = 1, size = 3) +
    scale_color_discrete(
      # #values = c("#7F3C8D",  "#11A579",  "#3969AC",  "#F2B701",  "#E73F74",  "#80BA5A",  "#E68310",  "#008695",  "#CF1C90",  "#F97B72",  "#4B4B8F",  "#A5AA99")
      # values = colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(length(unique(zi_pi$module)))
      )+#自定义颜色，按module为节点着色；本例只有9个modules,多提供了几种颜色示例。
    theme_bw()+
    theme(panel.grid = element_blank(), 
          axis.line = element_line(colour = 'black'),
          panel.background = element_blank(), 
          legend.key = element_blank(),
          legend.text = element_text(color="black",size=16,family = "serif",face = "bold"),#设置图例字体。serif表示是Times New Roman字体，想用Arial的话serif改sans，下文同。face = "bold"设置字体为粗体。如果不想粗体，bold改成plain
          legend.title = element_text(color="black",size=18,family = "serif",face = "bold"),#设置图例标题字体。
          axis.text = element_text(color="black",size=16,family = "serif",face = "bold"),#设置坐标轴刻度文本
          axis.title = element_text(color="black",size=18,family = "serif",face = "bold")) +#设置坐标轴标题字体。
    labs(x = 'Among-module connectivity (Pi)', y = 'Within-module connectivity (Zi)', color = 'Modules') +#设置x轴和y轴标题的名称，color = ''设置图例的标题名称
    #scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+#可自定义x轴范围及刻度
    #scale_y_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 1))+#可自定义y轴范围及刻度
    geom_vline(xintercept = 0.62, linetype="dotted",size = 0.6) +#在x轴0.62处添加竖直线
    geom_hline(yintercept = 2.5, linetype="dotted",size = 0.6) +#在y轴2.5处添加水平线
    if (lable) {
      geom_label_repel(data = zi_pi, aes(label = label),
                       size = 3,#标签字体大小
                       fontface="bold", #bold为粗体，plain为非粗体
                       box.padding = unit(0.5, "lines"),#标签文字到点的距离
                       point.padding = unit(0.8, "lines"),#点周围填充
                       segment.color = "black",#指示点的短线段的颜色
                       show.legend = FALSE, 
                       max.overlaps = 10000)
    }
}

#' ZIPI图绘制函数
#' @param zi_pi 标准ZIPI表格输入
#' @return 返回ggplot图形文件
zipi_draw <- function(zi_pi){
  library(ggrepel)#避免多标签时标签重叠
  library(readr)
  p2 <- ggplot(zi_pi, aes(p, z)) +
    # 注意，这里的roles如果为数字，aes上色会将其变为连续型变量，应该指定其为离散型
    geom_point(aes(color = module), alpha = 1, size = 3) +
    # scale_color_manual(values = c("#7F3C8D",  "#11A579",  "#3969AC",  "#F2B701",  "#E73F74",  "#80BA5A",  "#E68310",  "#008695",  "#CF1C90",  "#F97B72",  "#4B4B8F",  "#A5AA99"))+#自定义颜色，按module为节点着色；本例只有9个modules,多提供了几种颜色示例。
    scale_color_manual(
      values = colorRampPalette(RColorBrewer::brewer.pal(9,"Set2"))(length(unique(zi_pi$module)))
    )+ #自定义颜色
    theme_bw()+
    theme(panel.grid = element_blank(), 
          axis.line = element_line(colour = 'black'),
          panel.background = element_blank(), 
          legend.key = element_blank(),
          legend.text = element_text(color="black",size=12,family = "serif",face = "bold"),#设置图例字体。serif表示是Times New Roman字体，想用Arial的话serif改sans，下文同。face = "bold"设置字体为粗体。如果不想粗体，bold改成plain
          legend.title = element_text(color="black",size=14,family = "serif",face = "bold"),#设置图例标题字体。
          axis.text = element_text(color="black",size=12,family = "serif",face = "bold"),#设置坐标轴刻度文本
          axis.title = element_text(color="black",size=14,family = "serif",face = "bold")) +#设置坐标轴标题字体。
    labs(x = 'Among-module connectivity (Pi)', y = 'Within-module connectivity (Zi)', color = 'Modules') +#设置x轴和y轴标题的名称，color = ''设置图例的标题名称
    #scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+#可自定义x轴范围及刻度
    #scale_y_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 1))+#可自定义y轴范围及刻度
    geom_vline(xintercept = 0.62, linetype="dotted",size = 0.6) +#在x轴0.62处添加竖直线
    geom_hline(yintercept = 2.5, linetype="dotted",size = 0.6)#在y轴2.5处添加水平线
  p2
  return(p2)
}
# zipi_draw(zi_pi = zi_pi),taxa.roles))
}


