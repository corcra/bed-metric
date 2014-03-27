# Give this script a list of .metric outputs and it'll plot them all on the same graph!

args<-commandArgs(TRUE)
THRESHOLD<-0.01

library(ggplot2)
library(grid)

mytheme<-theme(legend.position="bottom",legend.direction="vertical",legend.key=element_rect(fill="white"),axis.title.x=element_text(vjust=-0.2),axis.title.y=element_text(angle=90,vjust=0.2),plot.title=element_text(size=15,vjust=1.4),plot.margin=unit(c(0.75,0.75,0.75,0.75),"cm"),panel.background=element_rect(fill="grey97",colour=NA),axis.line=element_line(colour="grey"))

get_deriv_model<-function(dataset){
    n<-nrow(dataset)
    deriv<-vector(length=n)
    prev_y<-0
    prev_x<-0
    before<-0
    after<-0
    crossed<-0
    x_val<-0
    y_val<-0
    model<-NA
    total_sorted<-sort(dataset[,4])
    uniq_sorted<-dataset[,3][order(dataset[,4])]

    for (i in seq(n)){
        y<-uniq_sorted[i]
        x<-total_sorted[i]
        deriv[i]<-(y-prev_y)/(x-prev_x)
        if ((deriv[i]<THRESHOLD)&(crossed==0)){
            crossed<-1
            before_x<-total_sorted[(i-1)]
            after_x<-total_sorted[i]
            before_y<-uniq_sorted[(i-1)]
            after_y<-uniq_sorted[i]
            x_val<-(0.5)*(after_x+before_x)
            y_val<-(0.5)*(after_y+before_y)
            }
        prev_y<-y
        prev_x<-x
        }

    if (x_val==0){
        model<-lm(deriv[(n-5):n]~total_sorted[(n-5):n])
        intercept<-summary(model)$coef[1]
        slope<-summary(model)$coef[2]
        x_val<-(THRESHOLD-intercept)/slope
        real_intercept<-max(uniq_sorted)-slope*max(total_sorted)^2*(0.5)-intercept*max(total_sorted)
        y_val<-(0.5)*slope*(x_val)^2+ intercept*x_val + real_intercept
        }
    return(list('deriv'=deriv,'x_val'=x_val,'y_val'=y_val,'lm'=model))
    }

colours<-rainbow(length(args),v=0.9)
colours_all<-c()
linetypes<-c()
i=0
total_all<-c()
uniq_all<-c()
name_all<-c()
name_order<-c()
asymptotes<-data.frame()
for (filename in args){
    cat(paste0("Loading ",filename,"!\n"))
    # how many data series have we included? (this is for getting colours)
    i = i+1
    # read data
    the_data<-read.table(filename,as.is=TRUE,skip=1)

    # name for the data (this is the filename after stripping directories etc., so you can give crazy long paths if you want)
    dataname<-gsub(".record","",strsplit(filename,'/')[[1]][length(strsplit(filename,'/')[[1]])])

    # these will be added to the total count, but first we have to check if they need to be artificially extended
    total<-the_data[,4]/1e6
    uniq<-the_data[,3]/1e6

    # true name
    name<-rep(dataname,nrow(the_data))
    name_order<-c(name_order,dataname)

    # add the colour
    colours_all<-c(colours_all,colours[i])

    # add linetype
    linetypes<-c(linetypes,"solid")

    # now, we check the derivative model to see if we need to predict any values
    deriv_data<-get_deriv_model(the_data)
    model<-deriv_data$lm

    # these will be added as dotted lines later
    x_metric<-deriv_data$x_val
    y_metric<-deriv_data$y_val

    if(!is.null(names(model))){
        # model is non-null, which means values WERE predicted!
        n_fake<-nrow(the_data) #totally arbitrary
        m<-max(total)*1e6
        mu<-max(uniq)*1e6
        x_fake<-seq(m/1e6,x_metric/1e6,length.out=n_fake)
        # quadratic fit
        a<-coef(summary(model))[2]
        b<-coef(summary(model))[1]
        c<-mu-(0.5)*a*m^2-b*m
        y_fake<-(0.5)*a*(1e6*x_fake)^2/1e6+b*x_fake+c/1e6
       
        # fake data gets a dashed line...
        linetypes<-c(linetypes,"dashed")

        # fake data gets a fake name...
        fake_name<-rep(paste0(dataname," (projected)"),n_fake)
        name_order<-c(name_order,fake_name)

        # repeat the colour
        colours_all<-c(colours_all,colours[i])

        # add this to what we had
        total<-c(total,x_fake)
        uniq<-c(uniq,y_fake)
        name<-c(name,fake_name)
    }

    asymptote_name<-paste0(dataname," (asymptote)")
    name<-c(name,asymptote_name)
    name_order<-c(name_order,asymptote_name)
    total<-c(total,0)
    uniq<-c(uniq,0)
    linetypes<-c(linetypes,"dotted")
    colours_all<-c(colours_all,colours[i])
    asymptotes<-rbind(asymptotes,c(x_metric/1e6,y_metric/1e6,i))
#    asymptotes[[length(asymptotes)+1]]<-list("x_metric"=x_metric,"y_metric"=y_metric,"colour"=colours[i])
#    asymptotes[[length(asymptotes)+1]]<-geom_segment(aes(x=0,xend=x_metric/1e6,y=y_metric/1e6,yend=y_metric/1e6,linetype=asymptote_name))
#    asymptotes<-c(asymptotes,geom_segment(aes(x=0,xend=x_metric/1e6,y=y_metric/1e6,yend=y_metric/1e6,linetype=asymptote_name)))

    total_all<-c(total_all,total)
    uniq_all<-c(uniq_all,uniq)
    name_all<-c(name_all,name)
}
names(asymptotes)<-c("x_metric","y_metric","col")

cat("Plotting!\n")
name_labels<-unique(name_order)
name_all<-factor(name_all,name_labels)
# what's going on here is we don't want to look at any of the asymptotes in the key
n_l<-c()
for (lab in name_labels){
    if(grepl("asymptote",lab)){
        n_l<-c(n_l,NULL)
    } else{
        n_l<-c(n_l,lab)
    }
}

d<-data.frame(total_all,uniq_all,name_all)
names(d)<-c("total","unique","dataset")

write.table(d,file="test.txt",col.names=T,row.names=F,quote=F)
ggplot(d,aes(x=total,y=unique,linetype=dataset,colour=dataset))+geom_line()+scale_color_manual(values=colours_all,breaks=n_l,labels=n_l)+guides(dataset=guide_legend())+mytheme+xlab("total (millions)")+ylab("unique (millions)")+scale_linetype_manual(values=linetypes,breaks=n_l,labels=n_l)+geom_segment(data=asymptotes,aes(x=0,xend=x_metric,y=y_metric,yend=y_metric,colour=levels(d$dataset)[which(grepl("asymptote",levels(d$dataset)))][col],linetype=levels(d$dataset)[which(grepl("asymptote",levels(d$dataset)))][col]))+scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand=c(0,0),limit=c(0,max(d$unique)*1.03))+ggtitle("Dataset quality comparison")

#+guides(dataset=guide_legend(title=NULL,byrow=FALSE),linetype=FALSE)
#+theme(legend.position = "bottom")+guide_legend(nrow=2)
ggsave("test.pdf")
