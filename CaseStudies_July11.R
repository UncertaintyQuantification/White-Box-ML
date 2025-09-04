library(RobustGaSP)
library(ggplot2)
library(patchwork)

# load('Datasets_Jun2.Rdata')

source('Phase_Functions.R')


# LOAD ALL DATASETS ###########
X_train <- as.matrix(read.csv('All_Datasets/X_train.csv'))
G_entries_train <- as.matrix(read.csv('All_Datasets/G_entries_train.csv'))

X_test_2d <- as.matrix(read.csv('All_Datasets/X_test_2d.csv'))
G_entries_test_2d <- as.matrix(read.csv('All_Datasets/G_entries_test_2d.csv'))

X_test <- as.matrix(read.csv('All_Datasets/X_test.csv'))
y_test <- as.vector(read.csv('All_Datasets/y_test.csv'))$x

X_test_grid <- read.csv('All_Datasets/X_test_grid.csv')
phase_record_grid <- read.csv('All_Datasets/phase_record_grid.csv')

X_train_paired <- read.csv('All_Datasets/X_train_paired.csv')
y_train_paired <- as.vector(read.csv('All_Datasets/y_train_paired.csv'))$x

X_train_nonpaired <- read.csv('All_Datasets/X_train_nonpaired.csv')
y_train_nonpaired <- as.vector(read.csv('All_Datasets/y_train_nonpaired.csv'))$x

all_X_train <- array(0, dim = c(500, 2, 10))
all_X_train_nonpaired <- array(0, dim = c(500, 13, 10))
all_G_entries_train <- array(0, dim = c(500, 600, 10))
all_phases <- matrix(0, nrow = 500, ncol = 10)

for(i in 1:10){
  filename_for_2d_training_data <- paste0('All_Datasets/all_X_train_',i,'.csv')
  filename_for_13d_training_data <- paste0('All_Datasets/all_X_train_nonpaired_',i,'.csv')
  filename_for_2d_training_data_response <- paste0('All_Datasets/all_G_entries_train_',i,'.csv')
  filename_for_13d_training_data_response <- paste0('All_Datasets/all_phases_',i,'.csv')
  
  all_X_train[,,i] <- as.matrix(read.csv(filename_for_2d_training_data))
  all_X_train_nonpaired[,,i] <- as.matrix(read.csv(filename_for_13d_training_data))
  all_G_entries_train[,,i] <- as.matrix(read.csv(filename_for_2d_training_data_response))
  all_phases[,i] <- as.vector(read.csv(filename_for_13d_training_data_response))$x
  
}


n_train <- 50
ppgp <- ppgasp(design = X_train[1:n_train,], response = log(G_entries_train[1:n_train,]),
               isotropic = F, nugget.est = T, method = 'mle', optimization = 'nelder-mead')

grid_size = 26

case_r_1 <- case_study_alpha_vary(ppgp, r = 0.5, grid_size = grid_size)
case_r_2 <- case_study_alpha_vary(ppgp, r = 2.0, grid_size = grid_size)
case_lb_1 <- case_study_alpha_vary(ppgp, lb = 1e-3, grid_size = grid_size)
case_lb_2 <- case_study_alpha_vary(ppgp, lb = 1, grid_size = grid_size)
case_chi_1 <- case_study_alpha_vary(ppgp, chiABN = 25.0, grid_size = grid_size)
case_chi_2 <- case_study_alpha_vary(ppgp, chiABN = 50.0, grid_size = grid_size)

plt_r1 <- plot_case_output(case_r_1, title = expression(r == 0.5))
plt_r2 <- plot_case_output(case_r_2, title = expression(r == 2.0))
plt_lb1 <- plot_case_output(case_lb_1, title = expression(l[B]/b == 10^{-3}))
plt_lb2 <- plot_case_output(case_lb_2, title = expression(l[B]/b == 1))
plt_chi1 <- plot_case_output(case_chi_1, title = expression(N[A] ~ chi[AB] == 25.0))
plt_chi2 <- plot_case_output(case_chi_2, title = expression(N[A] ~ chi[AB] == 50.0))

plt_r1 <- plt_r1 + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
plt_lb1 <- plt_lb1 + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                           axis.title.y = element_blank(), axis.text.y = element_blank())
plt_chi1 <- plt_chi1 + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                             axis.title.y = element_blank(), axis.text.y = element_blank())
plt_r2 <- plt_r2 #+ 
   #theme(#axis.ticks.x = element_line(linewidth=.5), 
                         #axis.text.x = element_text(angle=10, vjust = 1, hjust = 1))
plt_lb2 <- plt_lb2 + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                           #axis.ticks.x = element_line(linewidth=.5), 
                           #axis.text.x = element_text(angle=10, vjust = 1, hjust = 1)
                           )
plt_chi2 <- plt_chi2 + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                             #axis.ticks.x = element_line(linewidth=.5), 
                             #axis.text.x = element_text(angle=10, vjust = 1, hjust = 1)
                             )

#png('figures/cases.png', width = 1800, height = 1200)
pdf('figures/cases.pdf', width = 15, height = 10)
wrap_plots(plt_r1, plt_lb1, plt_chi1, 
           plt_r2, plt_lb2, plt_chi2, 
           ncol = 3, nrow = 2,
           guides = 'collect') &
  theme(legend.position = 'none',
        plot.subtitle = element_text(size = unit(22,'pt')),
        plot.title = element_text(size = unit(25,'pt'), margin = ggplot2::margin(t=0,r=0,l=0,b = 0, unit='pt')),
        axis.text = element_text(size = unit(22, 'pt')),
        axis.title = element_text(size=unit(25,'pt')),
        panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        #panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
        plot.margin = ggplot2::margin(t=.2,r=.6,b=.1,l=.6, unit='cm'))
dev.off()

pdf('figures/cases_no_lb.pdf', width = 10, height = 10)
wrap_plots(plt_r1, plt_chi1,
           plt_r2, plt_chi2,
           nrow = 2, ncol = 2,
           guides = 'collect') &
  theme(legend.position = 'none',
        plot.subtitle = element_text(size = unit(22,'pt')),
        plot.title = element_text(size = unit(25,'pt'), margin = ggplot2::margin(t=0,r=0,l=0,b = 0, unit='pt')),
        axis.text = element_text(size = unit(22, 'pt')),
        axis.title = element_text(size=unit(25,'pt')),
        panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        #panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
        plot.margin = ggplot2::margin(t=.2,r=.6,b=.1,l=.6, unit='cm'))
dev.off()

# bar plot

time_truth_case <- case_r_2$sim_time[[3]]
time_pred_case <- case_r_2$pred_time[[3]]

bar_plt_case_df <- data.frame(
  Method = c('Sim', 'Pred'),
  Time = c(time_truth_case, time_pred_case)
)

#png('figures/time_case_barplt.png', height = 600, width = 500)
pdf('figures/time_case_barplt.pdf', width = 4, height = 5)
ggplot(bar_plt_case_df, aes(x = Method, y = Time)) + 
  geom_bar(stat = 'identity', aes(fill = Method), width = .75) +
  scale_fill_manual(values = c("Pred" = "#2980B9", "Sim" = "#2980B9")) +
  geom_text(nudge_y = 12, label = round(bar_plt_case_df$Time, 2),
            size = 18, size.unit = 'pt') +
  scale_x_discrete(expand = c(.5,.5)) +
  scale_y_continuous(breaks = c(0,50,100,150,200,250,300,350)) +
  labs(y = 'Time (s)', title = 'Generation Time') +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = unit(18, 'pt')),
        axis.title.y = element_text(size = unit(20, 'pt')),
        plot.title = element_text(size = unit(20, 'pt')),
        # panel.grid.minor = element_blank(),
        # panel.grid.major.x = element_blank(),
        # panel.grid.major.y = element_line(linewidth = 1),
        panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
        legend.position = "none")
dev.off()

#save.image('CaseStudies_July11.Rdata')
#load('CaseStudies_July11.Rdata')
