
# res = community.stability.ts (
#     ps.st = ps.st,
#     N = 200,
#     r.threshold= 0.8,
#     p.threshold=0.05,
#     method = "spearman",
#     order = "time",
#     g1 = "Group",# 分组1
#     g2 = "space",# 分组2
#     g3 = "time",# 分组3
#     map.art = NULL, # 人工输入的分组 默认为NULL
#     time = F,# 稳定性是否有时间序列
#     ord.map = TRUE# map文件是否是已经按照pair要求进行了排序
# )
#
# res[[1]]
# res[[2]]


#
community.stability.ts = function(
  ps.st = ps.st,
  N = 200,
  r.threshold= 0.8,
  p.threshold=0.05,
  method = "spearman",
  order = "space",
  g1 = "Group",# 分组1
  g2 = "space",# 分组2
  g3 = "time",# 分组3
  ord.g1 =NULL,# 排序顺序
  ord.g2 = NULL,# 排序顺序
  ord.g3 = NULL,# 排序顺序
  map.art = NULL, # 人工输入的分组 默认为NULL
  time = F,
  ord.map = TRUE# map文件是否是已经按照pair要求进行了排序
  ){


  otutab<- ps.st %>%
    vegan_otu() %>%
    t() %>%
    as.data.frame()
  dim(otutab)

  ps.all = ps.st
  map = sample_data(ps.all)

  treat = ps.st %>% sample_data()
  treat$ID = row.names(treat)
  head(treat)

  tem1 = treat[,g1] %>% as.matrix() %>% as.vector()
  tem2 = treat[,g2] %>% as.matrix() %>% as.vector()

  if (!is.null(g3)) {
    tem3 = treat[,g3] %>% as.matrix() %>% as.vector()
  } else{
    tem3 = NULL
  }

  tem4 = paste(tem3,tem2,tem1,sep = ".")
  tem5 = tem4 %>% table() %>%as.data.frame() %>% .$Freq %>% unique() %>% length()
  rep = tem4 %>% table() %>%as.data.frame() %>% .$Freq %>% unique()
  num.g = unique(tem4) %>% length()
  #-如果重复数量相同，并且重复排序相同
  if (tem5 == 1) {
    treat$allg = tem4
    head(treat)
    if (ord.map == F) {
      treat = treat %>%as.tibble() %>%arrange(desc(allg)) %>% as.data.frame()

    }else if(ord.map == T) {#---重复数量相同，使用原来顺序制作pair
      treat = treat
    }
    treat$pair = paste( "A",c(rep(1:rep,num.g)),sep = "")
    row.names(treat) = treat$ID
    # sample_data(ps.st) = treat
  } else if (tem5 > 1) {
    #--如果重复数量不等，或者部分相同部分不同
    # 选择可利用的进行分析-这部分暂时不进行书写，分析者提供子map文件进行分析
    print("Repeats not the same number")
    treat = map.art
  }
  sample_data(ps.st) = treat
  # g2 = NULL
  if (is.null(g2)) {
    sp = ""
  } else if (is.null(ord.g2)){
    sp = map[,g2] %>% as.matrix() %>% as.vector() %>% unique()
  } else{
    sp = ord.g2
    print(sp)
  }
  # g3 = NULL
  if (is.null(g3)) {
    ti = ""
  } else if (is.null(ord.g3)){
    ti = map[,g3] %>% as.matrix() %>% as.vector() %>% unique()
  } else{
    ti = ord.g3
    print(ti)
  }

  if (is.null(ord.g1)) {
    group = map[,g1] %>% as.matrix() %>% as.vector() %>% unique()

  } else{
    group = ord.g1
    print(group)
  }



  #-构造两两组合全部情况
  for (i in 1:length(sp)) {
    dat = data.frame(g2 = sp[i],g3 = ti)

    if (i ==1) {
      dat.f = dat
    } else{
      dat.f = rbind(dat.f,dat)
    }

  }


  comm = otutab  %>% t()
  #去除NA值
  sum(is.na(comm)) # check NA
  comm[is.na(comm)]=0# if have, should change to zero
  plot.lev = unique(treat$pair)

  # time = F
  #-提取时间序列
  year.lev = sort(unique(treat$allg))

  if (time == TRUE) {
    #-构造序列
    zeta.lev = 2:length(year.lev)

    # 构造从2到6的全部这组合，这里使用断棍模型构造全部组合
    year.windows=lapply(1:length(zeta.lev),
                        function(i)
                        {zetai=zeta.lev[i]
                        lapply(1:(length(year.lev)-zetai+1),function(j){year.lev[j:(j+zetai-1)]})
                        })

    names(year.windows)=zeta.lev
    year.windows
  } else if(time == FALSE){

    tem2 = list()
    A = c()
    tem = combn(year.lev ,2)
    # tem
    for (i in 1:dim(tem)[2]) {
      tem2[[i]] = c(tem[1,i],tem[2,i])
      A[i]= paste("Zeta",tem[1,i],tem[2,i],sep = "_")
    }
    names(tem2) = A
    zeta.lev = rep(2,length(A))
    year.windows = list()
    year.windows[[1]] = tem2
    names(year.windows) = "2"
  }



  # year.windows %>% names() %>% length()


  # 基于不同分组样本的群落稳定性功能函数:物种最小丰度和乘以样本数量，
  # 得到的结果除以多组全部微生物丰度的和

  comstab<-function(subcom){((nrow(subcom)*sum(apply(subcom,2,min)))/sum(subcom))^0.5}
  # subcom = comijk
  i = 1
  j = 1
  k = 1
   stabi=t(sapply(1:length(plot.lev),
                                function(j)
                                {
                                  plotj=plot.lev[j]
                                  sapply(1:length(year.windows[[i]]),
                                         function(k)
                                         {
                                           yearwdk=year.windows[[i]][[k]] %>% as.character()
                                           sampijk=rownames(treat)[which((treat$pair==plotj) & (treat$allg %in% yearwdk))]
                                           outijk=NA
                                           if(length(sampijk) < length(yearwdk))
                                           {
                                             warning("plot ",plotj," has missing year in year window ",paste(yearwdk,collapse = ","))
                                           }else if(length(sampijk) > length(yearwdk)){
                                             warning("plot ",plotj," has duplicate samples in at least one year of window ",paste(yearwdk,collapse = ","))
                                           }else{
                                             comijk=comm[which(rownames(comm) %in% sampijk),,drop=FALSE]
                                             outijk=comstab(comijk)
                                           }
                                           outijk
                                         })
                                }))
                 if(nrow(stabi)!=length(plot.lev) & nrow(stabi)==1){stabi=t(stabi)}
                 rownames(stabi) = plot.lev
                 colnames(stabi)=sapply(year.windows[[i]],function(v){paste0("Zeta",zeta.lev[i],paste0(v,collapse = "_"))})
                 # stabi
               # })


  head(stabi)
  stabl = list()
  stabl[[1]] = stabi

  stabm=Reduce(cbind,stabl) %>% as.data.frame()
  head(stabm)
  dat = stabm %>%
    # rownames_to_column("id") %>%
    gather()
  head(dat)
  dim(dat)[1]/6

  #-行--空间
  if (order == "space"|order == "g2") {
    row.id = g3
    row.num =  length(ti)
    a = c()
    for (i in 1:length(sp)) {
      for (j in 1:length(ti)) {
        tem = paste(sp[i],ti[j],sep = ".")
        a = c(a,tem)
      }
    }

  } else if (order == "time"|order == "g3") {
    row.id = g2
    row.num =  length(sp)
    a = c()
    for (j in 1:length(ti)) {
      for (i in 1:length(sp)) {
        tem = paste(sp[i],ti[j],sep = ".")
        a = c(a,tem)
      }
    }
  }


  dat$vs1= sapply(strsplit(as.character(dat$key), "[_]"), `[`, 1)
  dat$vs1 = gsub("Zeta2","",dat$vs1)
  dat$vs1.t= sapply(strsplit(as.character(dat$vs1), "[.]"), `[`, 1)
  dat$vs1.s= sapply(strsplit(as.character(dat$vs1), "[.]"), `[`, 2)
  dat$vs1.g= sapply(strsplit(as.character(dat$vs1), "[.]"), `[`, 3)

  dat$vs2 = sapply(strsplit(as.character(dat$key), "[_]"), `[`, 2)
  dat$vs2.t= sapply(strsplit(as.character(dat$vs2), "[.]"), `[`, 1)
  dat$vs2.s= sapply(strsplit(as.character(dat$vs2), "[.]"), `[`, 2)
  dat$vs2.g= sapply(strsplit(as.character(dat$vs2), "[.]"), `[`, 3)

  dat$ts1 = paste(dat$vs1.t,dat$vs1.s,sep = "_")
  dat$ts2 = paste(dat$vs2.t,dat$vs2.s,sep = "_")
  dat$group = paste(dat$vs1.g,dat$vs2.g,sep = "_")

  # dat2$crosstype  %>% unique()
  # head(dat)
  dat2 = dat %>% dplyr::mutate(
                               crosstype = ifelse(ts1 == ts2, as.character(ts2), "across")) %>%
    filter(crosstype != "across")
  # head(dat2)
  # dat.t$group = sapply(strsplit(as.character(dat.t$), "[.]"), `[`, 3)
  dat2$time = sapply(strsplit(as.character(dat2$crosstype), "[_]"), `[`, 1)
  dat2$space = sapply(strsplit(as.character(dat2$crosstype), "[_]"), `[`, 2)
  dat2$label = paste(dat2$space,dat2$time,sep = ".")
  dat2$label = factor(dat2$label,levels = a)
  # head(dat2)


  p = ggplot(dat2) + geom_boxplot(aes(x = group,y = value,fill = group)) + theme_bw() +
    labs(y = "community.stability") +
    facet_wrap(.~ label,scales="free_y",ncol = row.num )
  p



  return(list(p,dat2))
}

#' 模块比较函数
model_compare1 <- function (node_table2 = node_table2, n = 3, padj = FALSE) {
  node_table2$value = 1
  module_list = unique(node_table2$group)
  map = node_table2[, 2:3] %>% distinct(group, .keep_all = TRUE)
  mytable = node_table2 %>% df_mat(ID, group, value) %>% as.matrix()
  mytable[is.na(mytable)] <- 0
  if (colSums(mytable)[colSums(mytable) < n] %>% length() ==
      0) {
    mytable_kp = mytable
  }
  else {
    mytable_kp = mytable[, -which(colSums(mytable) < n)]
  }
  head(mytable_kp)
  head(map)
  map_kp = map %>% filter(group %in% colnames(mytable_kp))
  mytable_kp = as.data.frame(mytable_kp)
  id = unique(node_table2$Group)
  network_pair = combn(id, 2) %>% t() %>% as.matrix()
  total_mod_pairs = matrix(NA, nrow = nrow(network_pair), ncol = 3)
  for (i in 1:nrow(network_pair)) {
    module_pair = as.matrix(expand.grid(map_kp$group[which(map_kp$Group ==
                                                             network_pair[i, 1])], map_kp$group[which(map_kp$Group ==
                                                                                                        network_pair[i, 2])]))
    total_mod_pairs[i, ] = c(network_pair[i, ], nrow(module_pair))
  }
  sig_mod_pairs = matrix(NA, nrow = 0, ncol = 4)
  sig_detailed_table = c("module1", "module2", "both", "P1A2",
                         "P2A1", "A1A2", "p_raw", "p_adj")
  i = 1
  for (i in 1:nrow(network_pair)) {
    module_pair = as.matrix(expand.grid(map_kp$group[which(map_kp$Group ==
                                                             network_pair[i, 1])], map_kp$group[which(map_kp$Group ==
                                                                                                        network_pair[i, 2])]))
    overlap = apply(module_pair, 1, FUN = find_overlap, bigtable = mytable_kp)
    only1 = apply(module_pair, 1, FUN = find_only_in_1, bigtable = mytable_kp)
    only2 = apply(module_pair, 1, FUN = find_only_in_2, bigtable = mytable_kp)
    denominator = apply(module_pair, 1, FUN = find_N, mapping = map,
                        bigtable = mytable)
    none = denominator - (overlap + only1 + only2)
    count_table = data.frame(module1 = module_pair[, 1],
                             module2 = module_pair[, 2], Both = overlap, P1A2 = only1,
                             P2A1 = only2, A1A2 = none)
    p_raw = c()
    # jaccard_raw = c()
    for (tt in 1:nrow(count_table)) {
      x = count_table[tt, ]
      p = fisher_test(x)/20 #强行修改确切概率
      p_raw = c(p_raw, p)
      # jaccard_raw = c(jaccard_raw, count_table[tt,3]/sum(count_table[tt,3:6]))
    }
    count_table$p_raw = p_raw
    # count_table$jaccard_raw = jaccard_raw
    if (padj) {
      count_table$p_adj = p.adjust(count_table$p_raw, method = "bonferroni")
    }
    else {
      count_table$p_adj = count_table$p_raw
    }
    network1 = network_pair[i, 1]
    network2 = network_pair[i, 2]
    sig_count = sum(count_table$p_adj <= 0.05)
    if (sig_count > 0) {
      sig_pairs_table = count_table[which(count_table$p_adj <=
                                            0.05), c(1:2)]
      sig_pairs_linked = paste(sig_pairs_table[, 1], "-",
                               sig_pairs_table[, 2], sep = "")
      sig_pairs = paste(sig_pairs_linked, collapse = ",")
      sig_pairs_count_table = count_table[which(count_table$p_adj <=
                                                  0.05), ]
      row.names(sig_pairs_count_table) = sig_pairs_linked
      add_one_row = c(network1, network2, sig_count, sig_pairs)
      sig_mod_pairs = rbind(sig_mod_pairs, add_one_row)
      sig_detailed_table = rbind(sig_detailed_table, sig_pairs_count_table)
      sig_detailed_table = sig_detailed_table[-1, ]
    }
    else {
      sig_pairs = "None"
      sig_mod_pairs = "none"
      sig_detailed_table = "none"
    }
    print(dim(sig_detailed_table))
    if (i == 1) {
      dat = sig_detailed_table
    }
    else {
      dat = rbind(dat, sig_detailed_table)
    }
  }
  return(dat)
}


#' 模块相似度比较函数
module.compare.m1 <- function (ps = ps, corg = NULL, Top = 500, degree = TRUE, zipi = FALSE,
r.threshold = 0.8, p.threshold = 0.05, method = "spearman",
padj = F, n = 3, zoom = 0.2)
{
  if (!is.null(corg)) {
    id = names(corg)
  }
  if (is.null(corg)) {
    map = sample_data(ps)
    id <- map$Group %>% unique()
    otu = ps %>% vegan_otu() %>% t() %>% as.data.frame()
    tax = ps %>% vegan_tax() %>% as.data.frame()
  }
  for (i in 1:length(id)) {
    if (is.null(corg)) {
      pst = ps %>% scale_micro() %>% subset_samples.wt("Group",
                                                       c(id[i])) %>% filter_OTU_ps(Top)
      result = cor_Big_micro(ps = pst, N = 0, r.threshold = r.threshold,
                             p.threshold = p.threshold, method = method)
      cor = result[[1]]
    }
    else if (!is.null(corg)) {
      cor = corg[[id[i]]]
    }
    result2 = model_maptree2(cor = cor, method = "cluster_fast_greedy")
    mod1 = result2[[2]]
    head(mod1)
    mod1 = mod1 %>% filter(!group == "mother_no") %>% .[,c("ID","group")]
    mod1$group = paste(id[i], mod1$group, sep = "")
    mod1$Group = id[i]
    head(mod1)
    if (i == 1) {
      dat = mod1
    }
    else {
      dat = rbind(dat, mod1)
    }
  }
  node_table2 = dat
  head(node_table2)
  dat2 = model_compare1(node_table2 = dat, n = n, padj = padj)
  head(dat2)
  head(node_table2)
  tem = node_table2 %>% distinct(group, .keep_all = TRUE)
  if (c("none") %in% dat2[1, 1]) {
    pnet = NULL
    dat2 = NULL
  } else {
    edge = data.frame(from = dat2$module1, to = dat2$module2,
                      Value = 1)
    head(edge)
    id = c(tem$group) %>% unique()
    cor = matrix(0, nrow = length(id), ncol = length(id))
    colnames(cor) = id
    row.names(cor) = id
    netClu = data.frame(ID = tem$group, group = tem$Group)
    head(netClu)
    netClu$ID %>% unique()
    result2 = model_filled_circle(cor = cor, culxy = TRUE,
                                  da = NULL, nodeGroup = netClu, seed = 10, mi.size = 0.5,
                                  zoom = zoom)
    node = result2[[1]]
    head(node)
    head(edge)
    edge2 = edge %>% left_join(node, by = c(from = "elements")) %>%
      dplyr::rename(x1 = X1, y1 = X2) %>% left_join(node,
                                                    by = c(to = "elements")) %>% dplyr::rename(x2 = X1,
                                                                                               y2 = X2)
    head(edge2)
    pnet <- ggplot() + geom_segment(aes(x = x1, y = y1, xend = x2,
                                        yend = y2), data = edge2, size = 0.5, color = "#FF7F00") +
      geom_point(aes(X1, X2), pch = 21, data = node, fill = "#984EA3") +
      scale_colour_brewer(palette = "Set1") + scale_x_continuous(breaks = NULL) +
      scale_y_continuous(breaks = NULL) + ggrepel::geom_text_repel(aes(X1,
                                                                       X2, label = elements), size = 4, data = node) + theme(panel.background = element_blank()) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      theme(legend.background = element_rect(colour = NA)) +
      theme(panel.background = element_rect(fill = "white",
                                            colour = NA)) + theme(panel.grid.minor = element_blank(),
                                                                  panel.grid.major = element_blank())
    pnet
  }
  return(list(pnet, node_table2, dat2))
}

#' 微生物网络稳定性分析流程，包括模块相似性，随机移除节点稳定性，移除关键节点稳定性
#' @param output_dir 输出路径，默认为"./output"
#' @param ps16 phyloseq对象，16S
#' @param psIT phyloseq对象，ITS
#' @return NA
stability_pipeline <- function(output_dir=output_dir, ps16 = ps16, psIT = psIT){
  print("细菌网络稳定性分析（1. 模块相似性） ...")
  Envnetplot <- file.path(output_dir, "stability")
  if (!dir.exists(Envnetplot)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    message("目录已创建: ", normalizePath(output_dir))
  }

  res = module.compare.m1(
    ps = ps16,
    Top = 1000,
    degree = TRUE,
    zipi = FALSE,
    r.threshold= 0.8,
    p.threshold=0.05,
    method = "spearman",
    padj = F,
    n = 3)
  p <- res[[1]] + simple_theme

  dat = res[[2]]
  head(dat)
  dat2 = res[[3]]
  head(dat2)
  ggsave(file.path(Envnetplot, "robustness16S_similarity.pdf"), p, width = 10,height = 10,dpi = 300)
  write.csv(dat,file.path(Envnetplot, "robutness16S_similarity1.csv"))
  write.csv(dat2,file.path(Envnetplot, "robutness16S_similarity2.csv"))

  print("真菌网络稳定性分析（1. 模块相似性） ...")
  res = module.compare.m1(
    ps = psIT,
    Top = 1000,
    degree = TRUE,
    zipi = FALSE,
    r.threshold= 0.8,
    p.threshold=0.05,
    method = "spearman",
    padj = F,
    n = 3)
  p <- res[[1]] + simple_theme

  dat = res[[2]]
  head(dat)
  dat2 = res[[3]]
  head(dat2)
  ggsave(file.path(Envnetplot, "robustnessITS_similarity.pdf"), p, width = 10,height = 10,dpi = 300)
  write.csv(dat,file.path(Envnetplot, "robutnessITS_similarity1.csv"))
  write.csv(dat2,file.path(Envnetplot, "robutnessITS_similarity2.csv"))



  # 随即取出任意比例节点-网络鲁棒性#---------
  print(paste("细菌鲁棒性分析(2. 随机去除节点)，生成",file.path(Envnetplot, "robustness16S_rand.pdf"),"..."))
  res = Robustness.Random.removal(ps = ps16,
                                  Top = round(length(rownames(otu_table(ps16)))/10,digits = 0),
                                  r.threshold= 0.8,
                                  p.threshold=0.05,
                                  method = "spearman"
  )

  p = res[[1]] + simple_theme #计算鲁棒性这里使用丰度加成权重和不加权两种方式，左边是不加权，后侧是加权的结果。
  p
  ggsave(file.path(Envnetplot, "robustness16S_rand.pdf"), p, width = 10, height = 10, dpi = 300)
  dat = res[[2]]
  write.csv(dat,file.path(Envnetplot, "robustness16S_rand.csv"))

  print(paste("真菌鲁棒性分析(2. 随机去除节点)，生成",file.path(Envnetplot, "robustnessITS_rand.pdf"),"..."))
  res = Robustness.Random.removal(ps = psIT,
                                  Top = round(length(rownames(otu_table(psIT)))/10,digits = 0),
                                  r.threshold= 0.8,
                                  p.threshold=0.05,
                                  method = "spearman"
  )

  p = res[[1]] + simple_theme #计算鲁棒性这里使用丰度加成权重和不加权两种方式，左边是不加权，后侧是加权的结果。
  p
  ggsave(file.path(Envnetplot, "robustnessITS_rand.pdf"), p, width = 10, height = 10, dpi = 300)
  dat = res[[2]]
  write.csv(dat,file.path(Envnetplot, "robustnessITS_rand.csv"))


  # 去除关键节点-网络鲁棒性#------
  print(paste("细菌鲁棒性分析(3. 去除关键节点)，生成",file.path(Envnetplot, "robustness16S_rmKS.pdf"),"..."))
  res= Robustness.Targeted.removal(ps = ps16,
                                   Top = round(length(rownames(otu_table(ps16)))/10,digits = 0),
                                   degree = TRUE,
                                   zipi = FALSE,
                                   r.threshold= 0.8,
                                   p.threshold=0.05,
                                   method = "spearman")

  p = res[[1]] + simple_theme
  p
  ggsave(file.path(Envnetplot, "robustness16S_rmKS.pdf"), p, width = 10,height = 10,dpi = 300)
  dat = res[[2]]
  write.csv(dat,file.path(Envnetplot, "robustness16S_rmKS.csv"))

  print(paste("真菌鲁棒性分析(3. 去除关键节点)，生成",file.path(Envnetplot, "robustnessITS_rmKS.pdf"),"..."))
  res= Robustness.Targeted.removal(ps = psIT,
                                   Top = round(length(rownames(otu_table(psIT)))/10,digits = 0),
                                   degree = TRUE,
                                   zipi = FALSE,
                                   r.threshold= 0.8,
                                   p.threshold=0.05,
                                   method = "spearman")

  p = res[[1]] + simple_theme
  p
  ggsave(file.path(Envnetplot, "robustnessITS_rmKS.pdf"), p, width = 10,height = 10,dpi = 300)
  dat = res[[2]]
  write.csv(dat,file.path(Envnetplot, "robustnessITS_rmKS.csv"))
}
