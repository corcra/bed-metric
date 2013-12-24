THRESHOLD<-0.01

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

the_data<-read.table(data,as.is=TRUE,skip=1)
deriv_data<-get_deriv_model(the_data)
deriv<-deriv_data$deriv
model<-deriv_data$lm
x_metric<-deriv_data$x_val
y_metric<-deriv_data$y_val
cat(paste(x_metric,y_metric,'\n'))

cat('Creating quality plot!\n')
lines(sort(the_data[,4])/1e6,the_data[order(the_data[,4]),3]/1e6,col=colour)
points(the_data[,4]/1e6,the_data[,3]/1e6,col=colour,cex=0.5,pch=20)
segments(x_metric/1e6,-5,x_metric/1e6,y_metric/1e6,lty=2,col=colour)
segments(-50,y_metric/1e6,x_metric/1e6,y_metric/1e6,lty=2,col=colour)
if(!is.null(names(model))){
    intercept<-max(the_data[,3])-coef(summary(model))[2]*max(the_data[,4])^2*(0.5)-coef(summary(model))[1]*max(the_data[,4])
    curve((coef(summary(model))[2]*(0.5)*(1e6*x)^2)/1e6+coef(summary(model))[1]*(x) + intercept/1e6,from=max(the_data[,4])/1e6,to=x_metric/1e6,add=T,lty=2,col=colour)
    }
