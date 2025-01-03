source( 'e1lib.R' ) 

parmmatrix
rownames(parmmatrix) <- c( 'λ0' , 'λ1' , 'µ0' , 'µ1' , 'q01' , 'q10'  )
pm <- parmmatrix 
pm
#       [,1]   [,2]   [,3]   [,4]   [,5]   [,6]   [,7]   [,8]   [,9]
# λ0  0.1800 0.1600 0.1400 0.1200 0.1000 0.0800 0.0600 0.0400 0.0200
# λ1  0.2000 0.2000 0.2000 0.2000 0.2000 0.2000 0.2000 0.2000 0.2000
# µ0  0.1000 0.1000 0.1000 0.1000 0.1000 0.1000 0.1000 0.1000 0.1000
# µ1  0.1000 0.1000 0.1000 0.1000 0.1000 0.1000 0.1000 0.1000 0.1000
# q01 0.0015 0.0015 0.0015 0.0015 0.0015 0.0015 0.0015 0.0015 0.0015
# q10 0.0030 0.0045 0.0060 0.0075 0.0090 0.0105 0.0120 0.0135 0.0150


fns <- paste0( 'e1-dtfb-',1:9,'.rds' )
alffns <- paste0( 'e1-dtfbalf-',1:9,'.rds' )

fbs <- lapply( fns, readRDS )
alffbs <- lapply( alffns, readRDS )


# function(fns)
for (i in seq_along(fbs)) fbs[[i]]$Experiment <- i 
for (i in seq_along(alffbs)) alffbs[[i]]$Experiment <- i 

do.call( rbind, lapply( fbs , function(fb )
# do.call( rbind, lapply( alffbs , function(fb)
{
	vnames <-  c('alpha', 'omega', 'yscale') 
	cx <- coef(fb)[vnames]
	odf = as.data.frame( cx )
	odf$lb <- exp( log(coef(fb)[vnames])-fb$err*1.96 )
	odf$ub <- exp( log(coef(fb)[vnames])+fb$err*1.96 )
	odf$Truth <- c( fb$theoralpha , fb$theoromega, NA )
	odf$Experiment <- fb$Experiment 
	odf$Variable <- c('alpha', 'omega', 'yscale' )
	odf$pa <- fb$theorpa 
	odf$empiricalpa <- fb$empiricalpa
	odf
	# fb$empiricalpa
	# fb$theorpa 
})) -> pldf 


# Filter to just alpha and omega rows and reshape data
pldf_alpha <- pldf[grepl("^alpha[0-9]*$", rownames(pldf)),]
pldf_omega <- pldf[grepl("^omega[0-9]*$", rownames(pldf)),]
# Combine into single dataframe
pldf1 <- data.frame(
  Experiment = pldf_alpha$Experiment,
  alpha_cx = pldf_alpha$cx,
  alpha_lb = pldf_alpha$lb,
  alpha_ub = pldf_alpha$ub,
  omega_cx = pldf_omega$cx,
  omega_lb = pldf_omega$lb,
  omega_ub = pldf_omega$ub,
  alpha_truth = pldf_alpha$Truth,
  omega_truth = pldf_omega$Truth
)
pcross <- ggplot(pldf1, aes(x = alpha_cx, y = omega_cx)) +
  # Error bars
  geom_errorbar(aes(ymin = omega_lb, ymax = omega_ub)) +
  geom_errorbarh(aes(xmin = alpha_lb, xmax = alpha_ub)) +
  # Estimated points
  geom_point(color = "blue", size = 2) +
  # Truth points
  geom_point(aes(x = alpha_truth, y = omega_truth), color = "red", size = 2) +
  # Faceting and formatting
  facet_wrap(~Experiment, nrow = 3) +
  theme_bw() +
  scale_x_continuous(limits = c(0, 12)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Alpha", y = "Omega")
pcross

pomega <- ggplot( pldf1, aes(x = omega_truth, y = omega_cx ) ) + 
	geom_errorbar( aes(ymin=omega_lb, ymax = omega_ub) , width=0) + 
	geom_point(size=4) + 
	geom_point(aes(x = omega_truth, y = omega_truth), colour = 'red', size=4) + 
	theme_bw() 
pomega


palpha <- ggplot( pldf1, aes(x = alpha_truth, y = alpha_cx ) ) + 
	geom_errorbar( aes(ymin=alpha_lb, ymax = alpha_ub), width=0 ) + 
	geom_point(size=4) + 
	geom_point(aes(x = alpha_truth, y = alpha_truth), colour = 'red', size=4) + 
	theme_bw() 
palpha

ggsave( pcross, file = 'e1plts-pcross.pdf' )
ggsave( pomega, file = 'e1plts-pomega.pdf' )
ggsave( pomega, file = 'e1plts-palpha.pdf' )

pldf
#                 cx         lb        ub Truth Experiment Variable
# alpha   1.77253487 1.25340671 2.5066723   2.0          1    alpha
# omega   0.71128857 0.61096530 0.8280854   0.9          1    omega
# yscale  0.86298620 0.80837494 0.9212868    NA          1   yscale
# alpha1  2.86286329 2.12830812 3.8509397   3.0          2    alpha
# omega1  0.68948518 0.59541474 0.7984179   0.8          2    omega
# yscale1 0.79132768 0.73929444 0.8470231    NA          2   yscale
# alpha2  3.37194285 2.68685622 4.2317108   4.0          3    alpha
# omega2  0.50198348 0.40401028 0.6237154   0.7          3    omega
# yscale2 0.81123686 0.75900719 0.8670606    NA          3   yscale
# alpha3  3.94082597 3.18998009 4.8684032   5.0          4    alpha
# omega3  0.46110103 0.36615335 0.5806697   0.6          4    omega
# yscale3 0.79258707 0.74099135 0.8477754    NA          4   yscale
# alpha4  4.86979293 4.03240324 5.8810793   6.0          5    alpha
# omega4  0.39969316 0.31317372 0.5101150   0.5          5    omega
# yscale4 0.75931858 0.70956605 0.8125596    NA          5   yscale
# alpha5  5.30136618 4.39639522 6.3926199   7.0          6    alpha
# omega5  0.28573570 0.20976074 0.3892287   0.4          6    omega
# yscale5 0.76574996 0.71602661 0.8189263    NA          6   yscale
# alpha6  5.22382865 4.29200554 6.3579568   8.0          7    alpha
# omega6  0.22880935 0.15876022 0.3297660   0.3          7    omega
# yscale6 0.77107797 0.72151842 0.8240417    NA          7   yscale
# alpha7  6.10423247 5.02195357 7.4197528   9.0          8    alpha
# omega7  0.15338727 0.09679879 0.2430573   0.2          8    omega
# yscale7 0.75458584 0.70554130 0.8070396    NA          8   yscale
# alpha8  6.09747706 4.93610345 7.5321003  10.0          9    alpha
# omega8  0.07668027 0.03804289 0.1545588   0.1          9    omega
# yscale8 0.77819584 0.72817724 0.8316502    NA          9   yscale
#                pa empiricalpa
# alpha   0.9148990       0.901
# omega   0.9148990       0.901
# yscale  0.9148990       0.901
# alpha1  0.8738980       0.873
# omega1  0.8738980       0.873
# yscale1 0.8738980       0.873
# alpha2  0.9062041       0.896
# omega2  0.9062041       0.896
# yscale2 0.9062041       0.896
# alpha3  0.8992824       0.891
# omega3  0.8992824       0.891
# yscale3 0.8992824       0.891
# alpha4  0.8894107       0.881
# omega4  0.8894107       0.881
# yscale4 0.8894107       0.881
# alpha5  0.8989883       0.890
# omega5  0.8989883       0.890
# yscale5 0.8989883       0.890
# alpha6  0.9075548       0.900
# omega6  0.9075548       0.900
# yscale6 0.9075548       0.900
# alpha7  0.9026192       0.898
# omega7  0.9026192       0.898
# yscale7 0.9026192       0.898
# alpha8  0.9106566       0.908
# omega8  0.9106566       0.908
# yscale8 0.9106566       0.908
