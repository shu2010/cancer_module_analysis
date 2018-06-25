library(shiny)
library(gridExtra)
library(ggnet)
library(ggplot2)
#library(sna)
#library(intergraph)
#library(network)
#load("/home/shu/VDLab_scripts/BioGrid/run_APcluster/Cluster_stage3/exmp_name.RData")

ui <- fluidPage(
  
  # App title ----
  titlePanel("Module consensus"),
  # selectInput("select", label = "Select Exemplar", 
  #             choices = c("select exemplar", c("TCGA_new_PRAD_3371", "TCGA_new_OV_2762", "TCGA_new_SARC_2120", 
  #                                              "TCGA_new_BRCA_985", "TCGA_new_OV_2951", "TCGA_new_PRAD_3255", 
  #                                              "TCGA_new_PRAD_3139", "TCGA_new_PAAD_3694", "TCGA_new_BRCA_543", 
  #                                              "TCGA_new_SKCM_1967", "TCGA_new_OV_2699", "TCGA_new_SARC_2128", 
  #                                              "TCGA_new_GBM_2492", "TCGA_new_LUAD_1047", "TCGA_new_LUAD_1417", 
  #                                              "TCGA_new_PRAD_3416", "TCGA_new_SARC_2214", "TCGA_new_SKCM_2002", 
  #                                              "TCGA_new_BRCA_691", "TCGA_new_GBM_2568", "TCGA_new_LUAD_1380", 
  #                                              "TCGA_new_PRAD_3320", "TCGA_new_LUAD_1109", "TCGA_new_SARC_2162", 
  #                                              "TCGA_new_BRCA_442", "TCGA_new_SARC_2255", "TCGA_new_BRCA_718", 
  #                                              "TCGA_new_BRCA_804", "TCGA_new_GBM_2476", "TCGA_new_SARC_2062", 
  #                                              "TCGA_new_PRAD_3167", "TCGA_new_PRAD_3312", "TCGA_new_PRAD_3487", 
  #                                              "TCGA_new_GBM_2555", "TCGA_new_PRAD_3575", "TCGA_new_SKCM_1647", 
  #                                              "TCGA_new_LUAD_1291", "TCGA_new_SARC_2192", "TCGA_new_BRCA_487", 
  #                                              "TCGA_new_GBM_2445", "TCGA_new_GBM_2388", "TCGA_new_GBM_2439", 
  #                                              "TCGA_new_SKCM_1736", "TCGA_new_SKCM_1757", "TCGA_new_SKCM_1874", 
  #                                              "TCGA_new_SKCM_1963", "TCGA_new_SKCM_1987", "TCGA_new_OV_2786", 
  #                                              "TCGA_new_PAAD_3754", "TCGA_new_SARC_2152", "TCGA_new_BRCA_416", 
  #                                              "TCGA_new_LUAD_1454", "TCGA_new_BRCA_278", "TCGA_new_LUAD_1423", 
  #                                              "TCGA_new_BRCA_293", "TCGA_new_BRCA_405", "TCGA_new_BRCA_932", 
  #                                              "TCGA_new_PRAD_3117", "TCGA_new_BRCA_650", "TCGA_new_BRCA_864", 
  #                                              "TCGA_new_LUAD_1132", "TCGA_new_PAAD_3601", "TCGA_new_PAAD_3731", 
  #                                              "TCGA_new_PRAD_3125", "TCGA_new_PRAD_3298", "TCGA_new_BRCA_172", 
  #                                              "TCGA_new_BRCA_558", "TCGA_new_BRCA_870", "TCGA_new_PRAD_3098", 
  #                                              "TCGA_new_PRAD_3392", "TCGA_new_PRAD_3446", "TCGA_new_SKCM_1808", 
  #                                              "TCGA_new_BRCA_164", "TCGA_new_BRCA_209", "TCGA_new_BRCA_813", 
  #                                              "TCGA_new_BRCA_688", "TCGA_new_LUAD_1021", "TCGA_new_LUAD_1240", 
  #                                              "TCGA_new_PAAD_3727", "TCGA_new_SARC_2163", "TCGA_new_GBM_2267", 
  #                                              "TCGA_new_GBM_2599", "TCGA_new_LUAD_1062", "TCGA_new_LUAD_1191", 
  #                                              "TCGA_new_OV_2796", "TCGA_new_PAAD_3750", "TCGA_new_BRCA_261", 
  #                                              "TCGA_new_BRCA_491", "TCGA_new_BRCA_716", "TCGA_new_PAAD_3746", 
  #                                              "TCGA_new_PRAD_3426", "TCGA_new_SKCM_1970", "TCGA_new_BRCA_871", 
  #                                              "TCGA_new_GBM_2426", "TCGA_new_OV_2916", "TCGA_new_PRAD_3369", 
  #                                              "TCGA_new_BRCA_977", "TCGA_new_LUAD_1178", "TCGA_new_PRAD_3191", 
  #                                              "TCGA_new_PRAD_3204", "TCGA_new_PRAD_3222", "TCGA_new_SARC_2051", 
  #                                              "TCGA_new_SKCM_1899", "TCGA_new_BRCA_704", "TCGA_new_BRCA_953", 
  #                                              "TCGA_new_OV_3046", "TCGA_new_PRAD_3367", "TCGA_new_SARC_2099", 
  #                                              "TCGA_new_SKCM_1928", "TCGA_new_BRCA_561", "TCGA_new_BRCA_608", 
  #                                              "TCGA_new_BRCA_896", "TCGA_new_PRAD_3085", "TCGA_new_PRAD_3228", 
  #                                              "TCGA_new_PRAD_3305", "TCGA_new_SKCM_1635", "TCGA_new_SKCM_1798", 
  #                                              "TCGA_new_BRCA_119", "TCGA_new_BRCA_304", "TCGA_new_BRCA_636", 
  #                                              "TCGA_new_BRCA_859", "TCGA_new_LUAD_1440", "TCGA_new_PRAD_3199", 
  #                                              "TCGA_new_PRAD_3473", "TCGA_new_BRCA_889", "TCGA_new_SARC_2158", 
  #                                              "TCGA_new_SKCM_1601", "TCGA_new_BRCA_344", "TCGA_new_BRCA_668", 
  #                                              "TCGA_new_PAAD_3724", "TCGA_new_PRAD_3364", "TCGA_new_BRCA_142", 
  #                                              "TCGA_new_BRCA_415", "TCGA_new_BRCA_837", "TCGA_new_PRAD_3114", 
  #                                              "TCGA_new_PRAD_3170", "TCGA_new_PRAD_3456", "TCGA_new_SARC_2022", 
  #                                              "TCGA_new_SKCM_1690", "TCGA_new_SKCM_1732", "TCGA_new_SARC_2030"
  #             )), 
  #             selected = NULL, multiple = FALSE),
  selectInput("select", label = "Select Exemplar", 
              choices = c("select exemplar", c("TCGA_new_BRCA_119", "TCGA_new_BRCA_142", "TCGA_new_BRCA_164", 
                                               "TCGA_new_BRCA_172", "TCGA_new_BRCA_209", "TCGA_new_BRCA_261", 
                                               "TCGA_new_BRCA_278", "TCGA_new_BRCA_293", "TCGA_new_BRCA_304", 
                                               "TCGA_new_BRCA_344", "TCGA_new_BRCA_405", "TCGA_new_BRCA_415", 
                                               "TCGA_new_BRCA_416", "TCGA_new_BRCA_442", "TCGA_new_BRCA_487", 
                                               "TCGA_new_BRCA_491", "TCGA_new_BRCA_543", "TCGA_new_BRCA_558", 
                                               "TCGA_new_BRCA_561", "TCGA_new_BRCA_608", "TCGA_new_BRCA_636", 
                                               "TCGA_new_BRCA_650", "TCGA_new_BRCA_668", "TCGA_new_BRCA_688", 
                                               "TCGA_new_BRCA_691", "TCGA_new_BRCA_704", "TCGA_new_BRCA_716", 
                                               "TCGA_new_BRCA_718", "TCGA_new_BRCA_804", "TCGA_new_BRCA_813", 
                                               "TCGA_new_BRCA_837", "TCGA_new_BRCA_859", "TCGA_new_BRCA_864", 
                                               "TCGA_new_BRCA_870", "TCGA_new_BRCA_871", "TCGA_new_BRCA_889", 
                                               "TCGA_new_BRCA_896", "TCGA_new_BRCA_932", "TCGA_new_BRCA_953", 
                                               "TCGA_new_BRCA_977", "TCGA_new_BRCA_985", "TCGA_new_GBM_2267", 
                                               "TCGA_new_GBM_2388", "TCGA_new_GBM_2426", "TCGA_new_GBM_2439", 
                                               "TCGA_new_GBM_2445", "TCGA_new_GBM_2476", "TCGA_new_GBM_2492", 
                                               "TCGA_new_GBM_2555", "TCGA_new_GBM_2568", "TCGA_new_GBM_2599", 
                                               "TCGA_new_LUAD_1021", "TCGA_new_LUAD_1047", "TCGA_new_LUAD_1062", 
                                               "TCGA_new_LUAD_1109", "TCGA_new_LUAD_1132", "TCGA_new_LUAD_1178", 
                                               "TCGA_new_LUAD_1191", "TCGA_new_LUAD_1240", "TCGA_new_LUAD_1291", 
                                               "TCGA_new_LUAD_1380", "TCGA_new_LUAD_1417", "TCGA_new_LUAD_1423", 
                                               "TCGA_new_LUAD_1440", "TCGA_new_LUAD_1454", "TCGA_new_OV_2699", 
                                               "TCGA_new_OV_2762", "TCGA_new_OV_2786", "TCGA_new_OV_2796", "TCGA_new_OV_2916", 
                                               "TCGA_new_OV_2951", "TCGA_new_OV_3046", "TCGA_new_PAAD_3601", 
                                               "TCGA_new_PAAD_3694", "TCGA_new_PAAD_3724", "TCGA_new_PAAD_3727", 
                                               "TCGA_new_PAAD_3731", "TCGA_new_PAAD_3746", "TCGA_new_PAAD_3750", 
                                               "TCGA_new_PAAD_3754", "TCGA_new_PRAD_3085", "TCGA_new_PRAD_3098", 
                                               "TCGA_new_PRAD_3114", "TCGA_new_PRAD_3117", "TCGA_new_PRAD_3125", 
                                               "TCGA_new_PRAD_3139", "TCGA_new_PRAD_3167", "TCGA_new_PRAD_3170", 
                                               "TCGA_new_PRAD_3191", "TCGA_new_PRAD_3199", "TCGA_new_PRAD_3204", 
                                               "TCGA_new_PRAD_3222", "TCGA_new_PRAD_3228", "TCGA_new_PRAD_3255", 
                                               "TCGA_new_PRAD_3298", "TCGA_new_PRAD_3305", "TCGA_new_PRAD_3312", 
                                               "TCGA_new_PRAD_3320", "TCGA_new_PRAD_3364", "TCGA_new_PRAD_3367", 
                                               "TCGA_new_PRAD_3369", "TCGA_new_PRAD_3371", "TCGA_new_PRAD_3392", 
                                               "TCGA_new_PRAD_3416", "TCGA_new_PRAD_3426", "TCGA_new_PRAD_3446", 
                                               "TCGA_new_PRAD_3456", "TCGA_new_PRAD_3473", "TCGA_new_PRAD_3487", 
                                               "TCGA_new_PRAD_3575", "TCGA_new_SARC_2022", "TCGA_new_SARC_2030", 
                                               "TCGA_new_SARC_2051", "TCGA_new_SARC_2062", "TCGA_new_SARC_2099", 
                                               "TCGA_new_SARC_2120", "TCGA_new_SARC_2128", "TCGA_new_SARC_2152", 
                                               "TCGA_new_SARC_2158", "TCGA_new_SARC_2162", "TCGA_new_SARC_2163", 
                                               "TCGA_new_SARC_2192", "TCGA_new_SARC_2214", "TCGA_new_SARC_2255", 
                                               "TCGA_new_SKCM_1601", "TCGA_new_SKCM_1635", "TCGA_new_SKCM_1647", 
                                               "TCGA_new_SKCM_1690", "TCGA_new_SKCM_1732", "TCGA_new_SKCM_1736", 
                                               "TCGA_new_SKCM_1757", "TCGA_new_SKCM_1798", "TCGA_new_SKCM_1808", 
                                               "TCGA_new_SKCM_1874", "TCGA_new_SKCM_1899", "TCGA_new_SKCM_1928", 
                                               "TCGA_new_SKCM_1963", "TCGA_new_SKCM_1967", "TCGA_new_SKCM_1970", 
                                               "TCGA_new_SKCM_1987", "TCGA_new_SKCM_2002")), 
  selected = NULL, multiple = FALSE),
  selectInput("Anno", label = "Annotation", choices = c("select Annotation", "MsigDB", "GO_MF", "GO_BP", "GO_CC"), 
              selected = NULL, multiple = FALSE),
  selectInput("Can", label = "Cancer", choices = c("Cancer type", c("BRCA", "GBM", "LUAD", "OV", "PAAD", "PRAD", "SARC", "SKCM", "All")), 
              selected = NULL, multiple = FALSE),
  # helpText("Degree based filtering, range: 100-50"),
  #   numericInput(inputId = "num", h3("Filtered modules"), n),
  sliderInput("num", "Consensus Module number", 
              min = 1, max = 7, value = 1, step= 1),
  sliderInput("in1", "Mutation frequency per patient", min = 1, max = 50, 
              value = 2, step= 1),
  
  #mainPanel(plotOutput(outputId = "distPlot", height = 1000, width = 1000))
  mainPanel(
    tabsetPanel(
      tabPanel("Plot",
               fluidRow(
                 plotOutput("distPlot"),
                 tableOutput("view")))),
    tabPanel("Summary",  tableOutput("view1"))
  )
)

#  )
##modules
#load("/home/shu/VDLab_scripts/BioGrid/run_APcluster/Cluster_stage3/trim_merge_modules_all.RData")
load("trim_merge_modules_all.RData")
#test_int <- lapply(trim_merge_modules_all, function(x)x[["TCGA_new_BRCA_261"]])[["TCGA_new_SARC_2120"]]
##plots
#load("/home/shu/VDLab_scripts/BioGrid/run_APcluster/Cluster_stage3/merge_vis_plots_all.RData")
load("merge_vis_plots_all.RData")
##Pathways
#load("/home/shu/VDLab_scripts/BioGrid/run_APcluster/Cluster_stage3/trim_merge_pathways_all.RData")
load("trim_merge_pathways_all.RData")
#d_tab <- trim_merge_pathways_max[[5]][["TCGA_new_SARC_2120"]]
#d_tab[,c(2,3)] <- sapply(d_tab[,c(2,3)], as.integer)
#samp_mut_mat
#load("~/VDLab_scripts/BioGrid/sample_mutation_matrix_union_biog.RData")
load("sample_mutation_matrix_union_biog_pub.RData")
#samp_mut_mat <- s3_mut_mat
#rm(s3_mut_mat)
##change sample name to current alias
##GO annotations
load("trim_merge_GO_MF_all.RData")
load("trim_merge_GO_BP_all.RData")
#trim_merge_GO_BP
load("trim_merge_GO_CC_all.RData")
##For cluster property plot
load("clust_pat_freq.RData")
##For survival plots
#load("survplots_merge_modules_all_alt.RData")
load("survplots_merge_modules_all_alt_new.RData")
surv_dat <- nt1
rm(nt1)

#t1 <- read.delim("~/VDLab_scripts/BioGrid/run_APcluster/cancer_id_pat_anno_alias_bg.tsv", header = T, sep= "\t")
#r1 <- as.character(t1[match(rownames(samp_mut_mat), t1$pat_id),3])
#rownames(samp_mut_mat) <- r1
#save(samp_mut_mat, file = "sample_mutation_matrix_union_biog_pub.RData")

# Define server logic ----
server <- function(input, output, session) { 
  
  # select <- reactive({
  #   validate(
  #     need(length(trim_merge_modules_all[[input$num]][[input$select]]) != 0, "Please select a different module")
  #   )
  #   get(trim_merge_modules_all[[input$num]][[input$select]])
  # })
  
  ##update number of proteins based on module number and cluster
  observeEvent(input$num,  {
    updateSliderInput(session = session, inputId = "in1", 
                      max = length(trim_merge_modules_all[[input$num]][[input$select]]))
  })
  
  output$distPlot <- renderPlot({
    validate(
      need(input$select != "select exemplar", "Please select an exemplar!!"),
      need(input$Anno != "select Annotation", "Please select annotation database!!"),
      need(input$Can != "Cancer type", "Please select Cancer type!!"),
      need(length(trim_merge_modules_all[[input$num]][[input$select]]) != 0, 
           "Please select a different module")
    )
    ##cancer frequency
    merge_modules_pat1 <- lapply(trim_merge_modules_all, function(x)x[[input$select]])
    # merge_modules_pat1 <- dataInput()
    merge_modules_pat2 <- lapply(merge_modules_pat1, function(x)rowSums(samp_mut_mat[,colnames(samp_mut_mat) %in% x]))
    
    can_freq <- lapply(merge_modules_pat2, function(x) x[x >= input$in1])
    
    can_freq1 <- lapply(can_freq,function(x)as.data.frame(table(unlist(lapply(strsplit(names(x), split = "_"), 
                                                                              function(y)y[3])))))
    # y <- select()
    # mut_can <- rowSums(samp_mut_mat[,colnames(samp_mut_mat) %in% y])
    # mut_can <- mut_can[mut_can >= input$in1 ]
    #  x <- as.data.frame(table(unlist(lapply(strsplit(names(mut_can), split = "_"), function(y)y[3]))))
    x <- can_freq1[[input$num]]
    #validate(need(!is.null(dim(x)), "Please select a lower frequency"))
    #need(!is.null(dim(x)), "Please select a lower frequency")  
    if(nrow(x) == 0){
      plot(1,1,col="white", axes=FALSE)
      text(1,1,"Please choose a lower frequency", cex = 2)
    }else{
      
      ##p_clust
      df_clust <- clust_pat_freq[[input$select]]
      colnames(df_clust) <- c("Cancer", "Patients")
      p_clu <- ggplot(data=df_clust, aes(x=Cancer, y=Patients)) + geom_bar(stat="identity", fill="gray") + 
        ggtitle(paste("Patients in Cluster", sum(df_clust$Patients), sep = ":")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text=element_text(size=8,face="bold"), plot.title = element_text(size = 10, face = "bold")) 
      
      #  annotate("text", x=dim(df_clust)[1]*0.9, y=max(df_clust$Patients)*0.9, 
      #          label= paste("Total Patients:", sum(df_clust$Patients), sep = "\n"), size = 3)
      p_clu <- p_clu + coord_flip() 
      
      ##p1 condition starts here
      colnames(x) <- c("Cancer", "Patients")
      
      #   plot(x, xlab = "Cancer", main = "Cancer prevalence")
      p <- ggplot(data=x, aes(x=Cancer, y=Patients)) + geom_bar(stat="identity", fill="gray") + 
        ggtitle(paste(input$in1, "Mutation(s) in", sum(x$Patients), "patients", sep = " ")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text=element_text(size=8,face="bold"), plot.title = element_text(size = 10, face = "bold")) 
      
      #  annotate("text", x=dim(x)[1]*0.9, y=max(x$Patients)*0.9, 
      #          label= paste("Total Patients:", sum(x$Patients), sep = "\n"), size = 3)
      p1 <- p + coord_flip() 
      
      ##mp1_prot
      
      
      #mat1 <- samp_mut_mat[,colnames(samp_mut_mat) %in% y]
      # mp1_prot <- colSums(mat1[rowSums(mat1) >= input$in1,,drop=FALSE])
      # d_mp1_prot <- as.data.frame(mp1_prot)
      trim_merge_modules_mat <- lapply(merge_modules_pat1, function(x)samp_mut_mat[,colnames(samp_mut_mat) %in% x])
      mp1_prot <- lapply(trim_merge_modules_mat, function(x) colSums(x[rowSums(x) >= input$in1,,drop=FALSE]))
      d_mp1_prot <- as.data.frame(mp1_prot[[input$num]])
      colnames(d_mp1_prot) <- "Frequency"
      d_mp1_prot$prot <- rownames(d_mp1_prot)
      p2 <- ggplot(data=d_mp1_prot, aes(x=prot, y=Frequency)) + geom_bar(stat="identity", fill="gray") + 
        ggtitle(paste("Patients:", sum(x$Patients), ";mut/pat:",input$in1, sep = "")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text=element_text(size=8,face="bold"), plot.title = element_text(size = 10, face = "bold")) 
      
      #    annotate("text", x=dim(d_mp1_prot)[1]*0.9, y=max(d_mp1_prot$Patients)*0.9, 
      #             label= paste("Total Patients:", sum(d_mp1_prot$Patients), sep = "\n"), size = 2.5)
      p2 <- p2 + coord_flip()
      ##total 
      
      #mat2_sum <- colSums(samp_mut_mat[,colnames(samp_mut_mat) %in% y])
      #d2 <- as.data.frame(mat2_sum)
      trim_merge_modules_freq <- lapply(merge_modules_pat1, function(x)colSums(samp_mut_mat[,colnames(samp_mut_mat) %in% x]))
      d2 <- as.data.frame(trim_merge_modules_freq[[input$num]])
      colnames(d2) <- "Frequency"
      d2$prot <- rownames(d2)
      
      pat_freq <- lapply(merge_modules_pat2, function(x) x[x >= 1])
      
      pat_freq1 <- lapply(pat_freq,function(x)as.data.frame(table(unlist(lapply(strsplit(names(x), split = "_"), 
                                                                                function(y)y[3])))))
      pat_freq2 <- pat_freq1[[input$num]]
      
      p3 <- ggplot(data=d2, aes(x=prot, y=Frequency)) + geom_bar(stat="identity", fill="gray") + 
        ggtitle(paste("Total patients", sum(pat_freq2$Freq), sep = ":")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text=element_text(size=8,face="bold"), plot.title = element_text(size = 10, face = "bold")) 
      
      #  annotate("text", x=dim(d2)[1]*0.9, y=max(d2$Patients)*0.9, 
      #          label= paste("Total Patients:", sum(d2$Patients), sep = "\n"), size = 2.5)
      p3 <- p3 + coord_flip() 
      ##graphs
      p4 <- merge_vis_plots_all[[input$num]][[input$select]]
      
      ##survival plots
      #Null plots
      if(is.null(surv_dat[[input$num]][[input$select]])){
      df <- data.frame()
      p5 <- ggplot(df) + geom_point() + xlim(0, 1) + ylim(0, 1) + ggtitle("No survival plot available for this module!!")
      }
      else if(is.null(surv_dat[[input$num]][[input$select]][[input$Can]])){
        df <- data.frame()
        p5 <- ggplot(df) + geom_point() + xlim(0, 1) + ylim(0, 1) + ggtitle("No survival plot available for this cancer type!!")
      }
      else if(!is.null(surv_dat[[input$num]][[input$select]][[input$Can]])){
        p5 <- surv_dat[[input$num]][[input$select]][[input$Can]]
      }
      grid.arrange(p_clu, p1, p2, p3, p4, p5, ncol=6, widths = c(5,5,5,5,8,8))
    }
  })
  ##table
  d_tab <- reactive({
    if(input$Anno == "MsigDB"){
      df_tab <- trim_merge_pathways_max[[input$num]][[input$select]]
      df_tab[,c(2,3)] <- sapply(df_tab[,c(2,3)], as.integer)
      df_tab
    }
    else if(input$Anno == "GO_MF"){
      df_tab <- trim_merge_GO_MF[[input$num]][[input$select]]
      df_tab
    }
    else if(input$Anno == "GO_BP"){
      df_tab <- trim_merge_GO_BP[[input$num]][[input$select]]
      df_tab
    }
    else if(input$Anno == "GO_CC"){
      df_tab <- trim_merge_GO_CC[[input$num]][[input$select]]
      df_tab
    }
  })
  output$view <- renderTable({
    d_tab()
  }, digits = 5)
  
}

# Run the app ----
shinyApp(ui = ui, server = server)