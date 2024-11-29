## 定义一个映射版本到可执行文件路径的列表
get_path <- function(version) {
  # 创建一个版本到路径的映射
  executable_paths <- list(
    "hiblup 1.5.3" = "./hiblup",  # 假设可执行文件名为hiblup
    "hiblup 1.5.3" = "./hiblup"   # 假设可执行文件名为hiblup
  )
  
  # 返回对应的可执行文件路径
  return(executable_paths[[version]])
}
HIBLUP_Test <-function(version_new, version_old){                                         ## 定义主函数  
  test_afreq <-function(version_new,version_old){                                         ## 定义等位基因频率的检测函数
    
    executable_new <- get_path(version_new)                                               ## 获得文件路径
    executable_old <- get_path(version_old)
    
    system(paste(executable_new, " --bfile demo/demo --allele-freq --out demo/demo_new"))  ## paste连接命令行 调用hiblup
    system(paste(executable_old, " --bfile demo/demo --allele-freq --out demo/demo_old"))
    Allele_new_freq <-read.table("demo/demo_new.afreq",header = TRUE)                      ## 读取文件
    Allele_old_freq <- read.table("demo/demo_old.afreq",header = TRUE)
     
    if (all.equal(Allele_new_freq, Allele_old_freq)) {                                     ## 判断语句 all.equal比较文件内容
      cat("\033[32m√\033[0m Normal Function")                                              ## 输出结果 
    } 
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
  }
  
  test_option_keep <- function(version_new,version_old){                                    ## 定义--keep参数的检测函数
    
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    
    system(paste(executable_new, " --bfile demo/demo --allele-freq --keep demo/id.filter.txt --out demo/id_new"))
    system(paste(executable_old, " --bfile demo/demo --allele-freq --keep demo/id.filter.txt --out demo/id_old"))
    ID_new_freq <-read.table("demo/id_new.afreq",header = TRUE) 
    ID_old_freq <- read.table("demo/id_old.afreq",header = TRUE)
    if (all.equal(ID_new_freq, ID_old_freq)) {
      cat("\033[32m√\033[0m Normal Function")
    } 
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
  }
  
  test_option_remove <- function(version_new,version_old){                      ## 定义参数 --remove检测函数
    
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    
    system(paste(executable_new," --bfile demo/demo --allele-freq --remove demo/id.filter.txt --out demo/idr_new"))
    system(paste(executable_old," --bfile demo/demo --allele-freq --remove demo/id.filter.txt --out demo/idr_old"))
    IDR_new_freq <-read.table("demo/idr_new.afreq",header = TRUE) 
    IDR_old_freq <- read.table("demo/idr_old.afreq",header = TRUE)
    if (all.equal(IDR_new_freq, IDR_old_freq)) {
      cat("\033[32m√\033[0m Normal Function")
    } 
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
  }
  
  test_option_extract <- function(version_new,version_old){                     ## 定义参数 --extract检测函数
    
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    
    system(paste(executable_new," --bfile demo/demo --allele-freq --extract demo/snp.filter.txt --out demo/snp_new"))
    system(paste(executable_old," --bfile demo/demo --allele-freq --extract demo/snp.filter.txt --out demo/snp_old"))
    SNP_new_freq <-read.table("demo/snp_new.afreq",header = TRUE)               
    SNP_old_freq <- read.table("demo/snp_old.afreq",header = TRUE)
    if (all.equal(SNP_new_freq, SNP_old_freq)) {
      cat("\033[32m√\033[0m Normal Function")
    } 
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
  }
  
  test_option_exclude <- function(version_new,version_old){                     ## 定义参数 --exclude检测函数
    
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    
    system(paste(executable_new," --bfile demo/demo --allele-freq --exclude demo/snp.filter.txt --out demo/snpr_new"))
    system(paste(executable_old," --bfile demo/demo --allele-freq --exclude demo/snp.filter.txt --out demo/snpr_old"))
    SNP_new_freq <-read.table("demo/snp_new.afreq",header = TRUE) 
    SNP_old_freq <- read.table("demo/snp_old.afreq",header = TRUE)
    if (all.equal(SNP_new_freq, SNP_old_freq)) {
      cat("\033[32m√\033[0m Normal Function")
    } 
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
  }
  
  test_gfreq <- function(version_new,version_old){                              ## 定义基因型频率检测函数
    
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    
    system(paste(executable_new," --bfile demo/demo --geno-freq --out demo/demo_new"))
    system(paste(executable_old," --bfile demo/demo --geno-freq --out demo/demo_old"))
    genotype_new_freq <-read.table("demo/demo_new.gfreq",header = TRUE)         
    genotype_old_freq <- read.table("demo/demo_old.gfreq",header = TRUE)
    if (all.equal(genotype_new_freq, genotype_old_freq)) {
      cat("\033[32m√\033[0m Normal Function")
    } 
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
  }
  
  test_transgeno<- function(version_new,version_old){                           ## 定义转换基因型格式检测函数
    
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    
    system(paste(executable_new," --trans-geno --bfile demo/demo --add --out demo/add_coding_new"))
    system(paste(executable_old," --trans-geno --bfile demo/demo --add --out demo/add_coding_old"))
    system(paste(executable_new," --trans-geno --bfile demo/demo --dom --out demo/dom_coding_new"))
    system(paste(executable_old," --trans-geno --bfile demo/demo --dom --out demo/dom_coding_old"))
    system(paste(executable_new," --trans-geno --bfile demo/demo --blupf90 --add --out demo/bf90_new"))
    system(paste(executable_old," --trans-geno --bfile demo/demo --blupf90 --add --out demo/bf90_old"))
    ADD_new <- read.table("demo/add_coding_new.geno.A.txt",header = TRUE)                    
    ADD_old <- read.table("demo/add_coding_old.geno.A.txt",header = TRUE)
    DOM_new <- read.table("demo/dom_coding_new.geno.D.txt",header = TRUE)
    DOM_old <- read.table("demo/dom_coding_new.geno.D.txt",header = TRUE)
    BF90_new <- read.table("demo/bf90_new.geno.A.txt",header = TRUE)
    BF90_old <- read.table("demo/bf90_old.geno.A.txt",header = TRUE)
    
    if (all.equal(ADD_new, ADD_old)) {
      cat("\033[32m√\033[0m Normal Function")
    } 
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    
    if (all.equal(DOM_new, DOM_old)) {
      cat("\033[32m√\033[0m Normal Function")
    } 
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    
    if (all.equal(BF90_new, BF90_old)) {
      cat("\033[32m√\033[0m Normal Function")
    } 
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
  }
  
  test_homo <- function(version_new,version_old){                               ## 定义纯合性检测函数
    
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    
    system(paste(executable_new," --bfile demo/demo --homo --out demo/demo_new"))
    system(paste(executable_old," --bfile demo/demo --homo --out demo/demo_old"))
    HOMO_new <- read.table("demo/demo_new.homo",header = TRUE)                     
    HOMO_old <- read.table("demo/demo_old.homo",header = TRUE)
    if (all.equal(HOMO_new, HOMO_old)) {
      cat("\033[32m√\033[0m Normal Function")
    } 
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }  
  }
  
  test_hete <- function(version_new,version_old){                               ## 定义杂合性检测函数
    
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    
    system(paste(executable_new," --bfile demo/demo --hete --out demo/demo_new"))
    system(paste(executable_old," --bfile demo/demo --hete --out demo/demo_old"))
    HETE_new <- read.table("demo/demo_new.hete",header = TRUE)                  
    HETE_old <- read.table("demo/demo_old.hete",header = TRUE)
    if (all.equal(HETE_new, HETE_old)) {
      cat("\033[32m√\033[0m Normal Function")
    } 
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }  
  }
  
  test_xrm <- function(version_new,version_old){                                ## 定义亲缘关系矩阵检测函数
    
    executable_new <-get_path(version_new)
    executable_old <-get_path(version_old)
    
    system(paste(executable_new," --make-xrm --pedigree demo/demo.ped --add --add-inv --write-txt --thread 32 --out demo/PRM_new"))
    system(paste(executable_old," --make-xrm --pedigree demo/demo.ped --add --add-inv --write-txt --thread 32 --out demo/PRM_old"))
    PA_new <- read.table("demo/PRM_new.txt")
    PA_old <- read.table("demo/PRM_old.txt")
    PAid_new <- read.table("demo/PRM_new.PA.id")
    PAid_old <- read.table("demo/PRM_old.PA.id")
    
    if (all.equal(PA_new,PA_old)){
      cat("\033[32m√\033[0m Normal Function")
    }
    else{
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    
    if (all.equal(PAid_new,PAid_old)){
      cat("\033[32m√\033[0m Normal Function")
    }
    else{
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    
    system(paste(executable_new," --make-xrm --bfile demo/demo --add  --add-inv --step 10000 --thread 32 --out demo/GRM_new --write-txt"))
    system(paste(executable_old," --make-xrm --bfile demo/demo --add  --add-inv --step 10000 --thread 32 --out demo/GRM_old --write-txt"))
    GRM_new <- read.table("demo/GRM_new.txt")
    GRM_old <- read.table("demo/GRM_old.txt")
    GAid_new <- read.table("demo/GRM_new.GA.id")
    GAid_old <- read.table("demo/GRM_old.GA.id")
    if (all.equal(GRM_new,GRM_old)){
      cat("\033[32m√\033[0m Normal Function")
    }
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    if (all.equal(GAid_new,GAid_old)){
      cat("\033[32m√\033[0m Normal Function")
    }
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    
    system(paste(executable_new," --make-xrm --bfile demo/demo --add --snp-weight demo/snp.weight.txt --thread 32 --out demo/wtGRM_new --write-txt"))
    system(paste(executable_old," --make-xrm --bfile demo/demo --add --snp-weight demo/snp.weight.txt --thread 32 --out demo/wtGRM_old --write-txt"))
    WTGRM_new <- read.table("demo/wtGRM_new.txt")
    WTGRM_old <- read.table("demo/wtGRM_old.txt")
    if (all.equal(WTGRM_new,WTGRM_old)){
      cat("\033[32m√\033[0m Normal Function")
    }
    else{
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    
    system(paste(executable_new," --make-xrm --pedigree demo/demo.ped --bfile demo/demo --alpha 0.05 --add --thread 32 --out demo/HRM_new --write-txt"))
    system(paste(executable_old," --make-xrm --pedigree demo/demo.ped --bfile demo/demo --alpha 0.05 --add --thread 32 --out demo/HRM_old --write-txt"))
    HRM_new <- read.table("demo/HRM_new.txt")
    HRM_old <- read.table("demo/HRM_old.txt")
    if (all.equal(HRM_new,HRM_old)){
      cat("\033[32m√\033[0m Normal Function")
    }
    else{
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    
    system(paste(executable_new," --make-xrm --pedigree demo/demo.ped --bfile demo/demo --alpha 0.05 --add-inv --thread 32 --out demo/Hinv_new --write-txt"))
    system(paste(executable_old," --make-xrm --pedigree demo/demo.ped --bfile demo/demo --alpha 0.05 --add-inv --thread 32 --out demo/Hinv_old --write-txt"))
    HI_new <- read.table("demo/Hinv_new.txt")
    HI_old <- read.table("demo/Hinv_old.txt")
    if (all.equal(HI_new,HI_old)){
      cat("\033[32m√\033[0m Normal Function")
    }
    else{
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    
    system(paste(executable_new," --make-xrm --pheno demo/demo.phe --rand 6,7 --out demo/Envrand_new --write-txt"))
    system(paste(executable_old," --make-xrm --pheno demo/demo.phe --rand 6,7 --out demo/Envrand_old --write-txt"))
    Envrand_new <- read.table("demo/Envrand_new.txt")
    Envrand_old<- read.table("demo/Envrand_old.txt")
    if (all.equal(Envrand_new,Envrand_old)){
      cat("\033[32m√\033[0m Normal Function")
    }
    else{
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    
  }
  
  test_ibc <- function(version_new,version_old){                                ## 定义近交系数检测函数
    
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    
    system(paste(executable_new," --ibc --pedigree demo/demo.ped --thread 32 --out demo/demo_new"))
    system(paste(executable_old," --ibc --pedigree demo/demo.ped --thread 32 --out demo/demo_old"))
    IBC_new <- read.table("demo/demo_new.ibc",header = TRUE)
    IBC_old <- read.table("demo/demo_old.ibc",header = TRUE)
    if (all.equal(IBC_new, IBC_old)) {
      cat("\033[32m√\033[0m Normal Function")
    } 
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    
    system(paste(executable_new," --ibc --bfile demo/demo  --thread 32 --out demo/GA_new"))
    system(paste(executable_old," --ibc --bfile demo/demo  --thread 32 --out demo/GA_old"))
    IBC1_new <- read.table("demo/GA_new.ibc",header = TRUE)
    IBC1_old <- read.table("demo/GA_old.ibc",header = TRUE)
    if (all.equal(IBC1_new, IBC1_old)) {
      cat("\033[32m√\033[0m Normal Function")
    } 
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
  }
  
  test_rc <- function(version_new,version_old){                                 ## 定义亲缘系数检测函数
    
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    
    system(paste(executable_new," --rc --pedigree demo/demo.ped --thread 32 --out demo/demo_new"))
    system(paste(executable_old," --rc --pedigree demo/demo.ped --thread 32 --out demo/demo_old"))
    RC_new <- read.table("demo/demo_new.rc",header = TRUE)
    RC_old <- read.table("demo/demo_old.rc",header = TRUE)
    if (all.equal(RC_new, RC_old)) {
      cat("\033[32m√\033[0m Normal Function")
    } 
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    
    system(paste(executable_new," --rc --bfile demo/demo  --thread 32 --out demo/GARC_new"))
    system(paste(executable_old," --rc --bfile demo/demo  --thread 32 --out demo/GARC_old"))
    RC1_new <- read.table("demo/GARC_new.rc",header = TRUE)
    RC1_old <- read.table("demo/GARC_old.rc",header = TRUE)
    if (all.equal(RC1_new, RC1_old)) {
      cat("\033[32m√\033[0m Normal Function")
    } 
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }  
  }
  
  test_pca <-function(version_new,version_old){                                 ## 定义主成分分析检测函数
    
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    
    system(paste(executable_new," --pca --bfile demo/demo  --npc 5 --out demo/demo_new"))
    system(paste(executable_old," --pca --bfile demo/demo  --npc 5 --out demo/demo_old"))
    
    PCA_new <- read.table("demo/demo_new.pc",header = TRUE)
    PCA_old <- read.table("demo/demo_old.pc",header = TRUE)
    if (all.equal(PCA_new, PCA_old)) {
      cat("\033[32m√\033[0m Normal Function")
    } 
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }  
    
    PCP_new <- read.table("demo/demo_new.pcp",header = TRUE,sep ="\t")
    PCP_old <- read.table("demo/demo_old.pcp",header = TRUE,sep ="\t")
    if (all.equal(PCP_new, PCP_old)) {
      cat("\033[32m√\033[0m Normal Function")
    } 
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }  
  }
  
  test_stm <-function(version_new,version_old){                                 ## 定义单性状模型检测函数
    
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    
    system(paste(executable_new," --make-xrm --pedigree demo/demo.ped --bfile demo/demo --alpha 0.05 --add  --dom --thread 32 --out demo/HRM_new --write-txt"))
    system(paste(executable_old," --make-xrm --pedigree demo/demo.ped --bfile demo/demo --alpha 0.05 --add  --dom --thread 32 --out demo/HRM_old --write-txt"))
    HRM_new <- read.table("demo/HRM_new.txt")
    HRM_old <- read.table("demo/HRM_old.txt")
    
    if(all.equal(HRM_new, HRM_old)) {
      
      system(paste(executable_new," --single-trait --pheno demo/demo.phe --pheno-pos 8 --dcovar 2,3 --qcovar 4,5 --rand 6,7 --xrm demo/HRM_new.HA,demo/HRM_new.HD --add --dom --threads 32 --out demo/stm_new"))
      system(paste(executable_old," --single-trait --pheno demo/demo.phe --pheno-pos 8 --dcovar 2,3 --qcovar 4,5 --rand 6,7 --xrm demo/HRM_new.HA,demo/HRM_new.HD --add --dom --threads 32 --out demo/stm_old"))
      VARS_new <- read.table("demo/stm_new.vars",header = TRUE)                 ## H矩阵 比较.vars .beta .anova .rand文件内容
      VARS_old <- read.table("demo/stm_old.vars",header = TRUE)                 
      BETA_new <- read.table("demo/stm_new.beta",header = TRUE)
      BETA_old <- read.table("demo/stm_old.beta",header = TRUE)
      ANOVA_new <- read.table("demo/stm_new.anova",header = TRUE)
      ANOVA_old <- read.table("demo/stm_old.anova",header = TRUE)
      RAND_new <- read.table("demo/stm_new.rand",header = TRUE)
      RAND_old <- read.table("demo/stm_old.rand",header = TRUE)
      if(all.equal(VARS_new, VARS_old)&all.equal(BETA_new, BETA_old)
         &all.equal(ANOVA_new, ANOVA_old)&all.equal(RAND_new, RAND_old)) {
        cat ("\033[32m√\033[0m Normal Function")
      }
      else{
        cat ("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
      }
    } 
    
    else {
      cat ("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    
  }
  
  test_rrstm <- function(version_new,version_old){                              ## 定义重复记录模型检测函数
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    
    system(paste(executable_new," --single-trait --pheno demo/demo.repeat.phe --pheno-pos 8 --rand 1 --bfile demo/demo --thread 32 --out demo/rr_new"))
    system(paste(executable_old," --single-trait --pheno demo/demo.repeat.phe --pheno-pos 8 --rand 1 --bfile demo/demo --thread 32 --out demo/rr_old"))
    VARS_new <- read.table("demo/rr_new.vars",header = TRUE)
    VARS_old <- read.table("demo/rr_old.vars",header = TRUE)
    BETA_new <- read.table("demo/rr_new.beta",header = TRUE)
    BETA_old <- read.table("demo/rr_old.beta",header = TRUE)
    RAND_new <- read.table("demo/rr_new.rand",header = TRUE)
    RAND_old <- read.table("demo/rr_old.rand",header = TRUE)
    if(all.equal(VARS_new, VARS_old)&all.equal(BETA_new, BETA_old)&all.equal(RAND_new, RAND_old)) {
      cat ("\033[32m√\033[0m Normal Function")
    }
    else{
      cat ("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
  }
  
  test_mtm <- function(version_new,version_old){                                ## 定义多性状模型  
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    
    system(paste(executable_new," --make-xrm --bfile demo/demo --add --dom --add-inv --step 10000 --thread 32 --out demo/GRM_new --write-txt"))
    system(paste(executable_old," --make-xrm --bfile demo/demo --add --dom --add-inv --step 10000 --thread 32 --out demo/GRM_old --write-txt"))
    GRM_new <- read.table("demo/GRM_new.txt")
    GRM_old <- read.table("demo/GRM_old.txt")
    GAid_new <- read.table("demo/GRM_new.GA.id")
    GAid_old <- read.table("demo/GRM_old.GA.id")
    if (all.equal(GRM_new,GRM_old)&all.equal(GAid_new,GAid_old)){
      system(paste(executable_new," --multi-trait --pheno demo/demo.phe --pheno-pos 8 9 10 --dcovar 2,3 0 2 --qcovar 4,5 5 4 --rand 6,7 7 0 --xrm demo/GRM_new.GA,demo/GRM_new.GD --vc-method AI --ai-maxit 30 --threads 32 --out demo/mtm_new"))
      system(paste(executable_old," --multi-trait --pheno demo/demo.phe --pheno-pos 8 9 10 --dcovar 2,3 0 2 --qcovar 4,5 5 4 --rand 6,7 7 0 --xrm demo/GRM_new.GA,demo/GRM_new.GD --vc-method AI --ai-maxit 30 --threads 32 --out demo/mtm_old"))
      VARS_new <- read.table("demo/mtm_new.vars",header = TRUE)                 ## 构建G矩阵 比较11对文件
      VARS_old <- read.table("demo/mtm_old.vars",header = TRUE)
      COVARS_new <- read.table("demo/mtm_new.covars",header = TRUE)
      COVARS_old <- read.table("demo/mtm_old.covars",header = TRUE)
      BETA1_new <- read.table("demo/mtm_new.T1.beta",header = TRUE)
      BETA1_old <- read.table("demo/mtm_old.T1.beta",header = TRUE)
      BETA2_new <- read.table("demo/mtm_new.T2.beta",header = TRUE)
      BETA2_old <- read.table("demo/mtm_old.T2.beta",header = TRUE)
      BETA3_new <- read.table("demo/mtm_new.T3.beta",header = TRUE)
      BETA3_old <- read.table("demo/mtm_old.T3.beta",header = TRUE)
      ANOVA1_new <- read.table("demo/mtm_new.T1.anova",header = TRUE)
      ANOVA1_old <- read.table("demo/mtm_old.T1.anova",header = TRUE)
      ANOVA2_new <- read.table("demo/mtm_new.T2.anova",header = TRUE)
      ANOVA2_old <- read.table("demo/mtm_old.T2.anova",header = TRUE)
      ANOVA3_new <- read.table("demo/mtm_new.T3.anova",header = TRUE)
      ANOVA3_old <- read.table("demo/mtm_old.T3.anova",header = TRUE)
      RAND1_new <- read.table("demo/mtm_new.T1.rand",header = TRUE)
      RAND1_old <- read.table("demo/mtm_old.T1.rand",header = TRUE)
      RAND2_new <- read.table("demo/mtm_new.T2.rand",header = TRUE)
      RAND2_old <- read.table("demo/mtm_old.T2.rand",header = TRUE)
      RAND3_new <- read.table("demo/mtm_new.T3.rand",header = TRUE)
      RAND3_old <- read.table("demo/mtm_old.T3.rand",header = TRUE)
      if(all.equal(VARS_new, VARS_old)&all.equal(COVARS_new, COVARS_old)&all.equal(BETA1_new, BETA1_old)
         &all.equal(ANOVA1_new, ANOVA1_old)&all.equal(RAND1_new, RAND1_old)&all.equal(BETA2_new, BETA2_old)
         &all.equal(ANOVA2_new, ANOVA2_old)&all.equal(RAND2_new, RAND2_old)&all.equal(BETA3_new, BETA3_old)
         &all.equal(ANOVA3_new, ANOVA3_old)&all.equal(RAND3_new, RAND3_old)) {
        cat ("\033[32m√\033[0m Normal Function")
      }
      else{
        cat ("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
      }
      
    }
    else {
      cat("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    
  }
  
  test_mme <- function(version_new,version_old){                                ## 定义混合线性模型求解的检测函数
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    
    system(paste(executable_new," --mme --pheno demo/demo.phe --pheno-pos 8 --rand 6 --bfile demo/demo --add --vc-priors 1,2,3 --pcg --threads 32 --out demo/mme_new"))
    system(paste(executable_old," --mme --pheno demo/demo.phe --pheno-pos 8 --rand 6 --bfile demo/demo --add --vc-priors 1,2,3 --pcg --threads 32 --out demo/mme_old"))
    RAND_new <- read.table("demo/mme_new.rand",header = TRUE)
    RAND_old <- read.table("demo/mme_old.rand",header = TRUE)
    BETA_new <- read.table("demo/mme_new.beta",header = TRUE)
    BETA_old <- read.table("demo/mme_old.beta",header = TRUE)
    if(all.equal(RAND_new, RAND_old)&all.equal(BETA_new, BETA_old)) {
      cat ("\033[32m√\033[0m Normal Function")
    }
    else{
      cat ("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
  }
  
  test_stat <- function(version_new,version_old){                               ## 定义随机效应显著性检验
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    
    system(paste(executable_new," --single-trait --pheno demo/demo.phe --pheno-pos 8 --rand 6,7 --bfile demo/demo --add --threads 32 --out demo/stat_new"))
    system(paste(executable_old," --single-trait --pheno demo/demo.phe --pheno-pos 8 --rand 6,7 --bfile demo/demo --add --threads 32 --out demo/stat_old"))
    system(paste(executable_new," --single-trait --pheno demo/demo.phe --pheno-pos 8 --bfile demo/demo --add --threads 32 --out demo/statt_new"))
    system(paste(executable_old," --single-trait --pheno demo/demo.phe --pheno-pos 8 --bfile demo/demo --add --threads 32 --out demo/statt_old"))
    stat_content <- readLines("demo/stat_new.log")                              ## 按行读取.log文件
    pattern <- "^\\[AI\\]"                                                              
    matches <- grep(pattern, stat_content)                                      ## 获得 文件中 含有 AI 行的索引
    ai_lines <- stat_content[matches]                                           ## 提取 文件中含有AI的行
    parsed_data <- lapply(ai_lines, function(line) {                            ## 提取并返回数值向量
      parts <- unlist(strsplit(line, "\\s+"))
      return(as.numeric(parts))
    })
    column_names <- c("Alg.","Iter.", "LogL.", "V(location)", "V(dam)","V(GA)","V(e)")
    data_frame <- do.call(rbind, parsed_data)                                   ## 合并数据
    colnames(data_frame) <- column_names
    demo <- as.data.frame(data_frame)
    
    statt_content <- readLines("demo/statt_new.log")
    pattern <- "^\\[AI\\]"
    matches <- grep(pattern, statt_content)
    ai_lines <- statt_content[matches]
    data <- lapply(ai_lines, function(line) {
      parts <- unlist(strsplit(line, "\\s+"))
      return(as.numeric(parts))
    })
    column_names <- c("Alg.","Iter.", "LogL.", "V(location)", "V(dam)")
    data_frame <- do.call(rbind, data)
    colnames(data_frame) <- column_names
    stat <- as.data.frame(data_frame)
    
    L1<-demo[nrow(demo),3]
    L0<-stat[nrow(stat),3]
    stat <- 2*(L1-L0)
    pvaue_new <- pchisq(stat,df = 1,lower.tail = FALSE)
    print(pvaue_new)
    
    stat_content <- readLines("demo/stat_old.log")
    pattern <- "^\\[AI\\]"
    matches <- grep(pattern, stat_content)
    ai_lines <- stat_content[matches]
    parsed_data <- lapply(ai_lines, function(line) {
      parts <- unlist(strsplit(line, "\\s+"))
      return(as.numeric(parts))
    })
    column_names <- c("Alg.","Iter.", "LogL.", "V(location)", "V(dam)","V(GA)","V(e)")
    data_frame <- do.call(rbind, parsed_data)
    colnames(data_frame) <- column_names
    demo <- as.data.frame(data_frame)
    
    statt_content <- readLines("demo/statt_old.log")
    pattern <- "^\\[AI\\]"
    matches <- grep(pattern, statt_content)
    ai_lines <- statt_content[matches]
    data <- lapply(ai_lines, function(line) {
      parts <- unlist(strsplit(line, "\\s+"))
      return(as.numeric(parts))
    })
    column_names <- c("Alg.","Iter.", "LogL.", "V(location)", "V(dam)")
    data_frame <- do.call(rbind, data)
    colnames(data_frame) <- column_names
    stat <- as.data.frame(data_frame)
    
    L1<-demo[nrow(demo),3]
    L0<-stat[nrow(stat),3]
    stat <- 2*(L1-L0)
    pvaue_old <- pchisq(stat,df = 1,lower.tail = FALSE)
    print(pvaue_old)
    
    if (pvaue_new == pvaue_old){
      cat ("\033[32m√\033[0m Normal Function")
    }
    else {
      cat ("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    
  }
  
  test_r2 <- function(version_new,version_old){                                 ## 定义可靠性检测函数
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    system(paste(executable_new," --single-trait --pheno demo/demo.phe --pheno-pos 8 --dcovar 2,3 --qcovar 4,5 --rand 7 --bfile demo/demo --r2 --threads 32 --out demo/r2_new"))
    system(paste(executable_old," --single-trait --pheno demo/demo.phe --pheno-pos 8 --dcovar 2,3 --qcovar 4,5 --rand 7 --bfile demo/demo --r2 --threads 32 --out demo/r2_old"))
    R2_new <- read.table("demo/r2_new.rand",header = TRUE)
    R2_old <- read.table("demo/r2_old.rand",header = TRUE)
    if (isTRUE(all.equal(R2_new,R2_old))){
      cat ("\033[32m√\033[0m Normal Function")
    }
    else {
      cat ("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
  }
  
  test_snpeff <- function(version_new,version_old){                             ## 定义标记效应检测函数
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    system(paste(executable_new," --single-trait --pheno demo/demo.phe --pheno-pos 8 --bfile demo/demo --pedigree demo/demo.ped --add --snp-effect --threads 32 --out demo/snp_new"))
    system(paste(executable_old," --single-trait --pheno demo/demo.phe --pheno-pos 8 --bfile demo/demo --pedigree demo/demo.ped --add --snp-effect --threads 32 --out demo/snp_old"))
    SNP_new <- read.table("demo/snp_new.snpeff",header = TRUE) 
    SNP_old <- read.table("demo/snp_old.snpeff",header = TRUE) 
    if (all.equal(SNP_new,SNP_old)){
      cat ("\033[32m√\033[0m Normal Function")
    }
    else{
      cat ("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
  }
  
  test_gebv <- function(version_new,version_old){                               ## 定义基因组育种值检测函数
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    system(paste(executable_new," --single-trait --pheno demo/demo.phe --pheno-pos 8 --bfile demo/demo --pedigree demo/demo.ped --add --snp-effect --threads 32 --out demo/snp_new"))
    system(paste(executable_old," --single-trait --pheno demo/demo.phe --pheno-pos 8 --bfile demo/demo --pedigree demo/demo.ped --add --snp-effect --threads 32 --out demo/snp_old"))
    SNP_new <- read.table("demo/snp_new.snpeff",header = TRUE) 
    SNP_old <- read.table("demo/snp_old.snpeff",header = TRUE) 
    if (isTRUE(all.equal(SNP_new,SNP_old))){
      system(paste(executable_new," --pred --bfile demo/demo --score demo/snp_new.snpeff --thread 10 --out demo/gebv_new"))
      system(paste(executable_old," --pred --bfile demo/demo --score demo/snp_old.snpeff --thread 10 --out demo/gebv_old"))
      GEBV_new <- read.table("demo/gebv_new.bv",header = TRUE)
      GEBV_old <- read.table("demo/gebv_old.bv",header = TRUE)
      if (all.equal(GEBV_new,GEBV_old)){
        cat ("\033[32m√\033[0m Normal Function")
        
      }
      else {
        cat ("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
      }
    }
    else{
      cat ("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    
    
  }
  
  test_gmat <- function(version_new,version_old){                               ## 定义基因组选配模块检测函数
    
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    system(paste(executable_new," --single-trait --pheno demo/demo.phe --pheno-pos 8 --bfile demo/demo --pedigree demo/demo.ped --add --snp-effect --threads 32 --out demo/snp_new"))
    system(paste(executable_old," --single-trait --pheno demo/demo.phe --pheno-pos 8 --bfile demo/demo --pedigree demo/demo.ped --add --snp-effect --threads 32 --out demo/snp_old"))
    SNP_new <- read.table("demo/snp_new.snpeff",header = TRUE) 
    SNP_old <- read.table("demo/snp_old.snpeff",header = TRUE) 
    if (all.equal(SNP_new,SNP_old)){
      system(paste(executable_new," --mating --bfile demo/demo --score demo/snp_new.snpeff --thread 32 --out demo/gmat_new"))
      system(paste(executable_old," --mating --bfile demo/demo --score demo/snp_old.snpeff --thread 32 --out demo/gmat_old"))
      GMAT_new <- read.table("demo/gmat_new.mating",header = TRUE)
      GMAT_old <- read.table("demo/gmat_old.mating",header = TRUE)
      if (all.equal(GMAT_new,GMAT_old)){
        cat ("\033[32m√\033[0m Normal Function")
      }
      else{
        cat ("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
      }
    }
    else{
      cat ("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
  }
  
  test_ld <- function(version_new,version_old){                                 ## 定义ld模块检测函数
    
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    system(paste(executable_new," --ld --bfile demo/demo --window-bp 1000000 --threads 32 --out demo/ld_new"))
    system(paste(executable_old," --ld --bfile demo/demo --window-bp 1000000 --threads 32 --out demo/ld_old"))
    LD_new <- read.table("demo/ld_new.info",header = TRUE)
    LD_old <- read.table("demo/ld_old.info",header = TRUE)
    if (all.equal(LD_new, LD_old)){
      cat ("\033[32m√\033[0m Normal Function")
    }
    else{
      cat ("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    system(paste(executable_new," --ld --bfile demo/demo --window-num 500 --threads 32 --out demo/ld1_new"))
    system(paste(executable_old," --ld --bfile demo/demo --window-num 500 --threads 32 --out demo/ld1_old"))
    LD1_new <- read.table("demo/ld1_new.info",header = TRUE)
    LD1_old <- read.table("demo/ld1_old.info",header = TRUE)
    if (all.equal(LD1_new, LD1_old)){
      cat ("\033[32m√\033[0m Normal Function")
    }
    else{
      cat ("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
    system(paste(executable_new," --ld --bfile demo/demo --window-geno --threads 32 --out demo/ld2_new"))
    system(paste(executable_old," --ld --bfile demo/demo --window-geno --threads 32 --out demo/ld2_old"))
    LD2_new <- read.table("demo/ld2_new.info",header = TRUE)
    LD2_old <- read.table("demo/ld2_old.info",header = TRUE)
    if (all.equal(LD_new, LD_old)){
      cat ("\033[32m√\033[0m Normal Function")
    }
    else{
      cat ("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
  }
  
  test_ldsc <- function(version_new,version_old){                               ## 定义ldsc模块检测函数
    executable_new <- get_path(version_new)
    executable_old <- get_path(version_old)
    system(paste(executable_new," --ldscore --bfile demo/demo --window-bp 1000000 --threads 32 --out demo/ldsc_new"))
    system(paste(executable_old," --ldscore --bfile demo/demo --window-bp 1000000 --threads 32 --out demo/ldsc_old"))
    LDSC_new <- read.table("demo/ldsc_new.ldsc",header = TRUE)
    LDSC_old <- read.table("demo/ldsc_old.ldsc",header = TRUE)
    if (all.equal(LDSC_new, LDSC_old)){
      cat ("\033[32m√\033[0m Normal Function")
    }
    else{
      cat ("The elements in the two files are not completely identical.please check your command or try to connect hibulp@author：Professor Lilin Yin.\n")
    }
  }
  ## 调用子检测函数
  test_afreq(version_new, version_old)
  test_option_keep(version_new, version_old)
  test_option_remove(version_new, version_old)
  test_option_extract(version_new, version_old)
  test_option_exclude(version_new, version_old)
  test_gfreq(version_new, version_old)
  test_transgeno(version_new, version_old)
  test_homo(version_new, version_old)
  test_hete(version_new, version_old)
  test_xrm(version_new, version_old)
  test_ibc(version_new, version_old)
  test_rc(version_new, version_old)
  test_pca(version_new, version_old)
  test_stm(version_new, version_old)
  test_rrstm(version_new, version_old)
  test_mtm(version_new, version_old)
  test_mme(version_new, version_old)
  test_stat(version_new, version_old)
  test_r2(version_new, version_old)
  test_snpeff(version_new, version_old)
  test_gebv(version_new, version_old)
  test_gmat(version_new, version_old)
  test_ld(version_new, version_old)
  test_ldsc(version_new, version_old)
}
## 调用主函数 检测功能是否正常
HIBLUP_Test("hiblup 1.5.3","hiblup 1.5.3")
system.time({HIBLUP_Test("hiblup 1.5.3","hiblup 1.5.3")})