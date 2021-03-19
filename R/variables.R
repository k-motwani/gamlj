Variables <- R6::R6Class("variables",
                            public=list(
                                dep=list(),
                                covs=list(),
                                factors=list(),
                                clusters=list(),
                                vars=list(),
                                options=list(),
                                data64=list(),
                                names64= names64$new(),
                                initialize = function(options) {
                                            if("cluster" %in% names(options)) clusters<-options$cluster else clusters<-NULL
                                             private$.initOptions(options)
                                             private$.init()
                                             
                                          },
                                fillinfo=function(data) {
                                  self$data64<-as.data.frame(data[0,self$vars])
                                  
                                  names(self$data64)<-jmvcore::toB64(self$vars)
                                  
                                  for (.var in self$factors) {
                                    if (!is.factor(data[[.var$name]])) data[[.var$name]]<-factor(data[[.var$name]])
                                    .var$levels<-levels(data[[.var$name]])
                                    .var$ncontrasts<-length(.var$levels)-1
                                    .var$cnames64<-unlist(lapply(seq_len(.var$ncontrast),function(i) paste0(.var$name64,"_._._",i)))
                                    .var$cnames<-unlist(lapply(seq_len(.var$ncontrast),function(i) paste0(.var$name,i)))
                                     stats::contrasts(self$data64[[.var$name64]]) <- lf.createContrasts(.var$levels, .var$contrast)
                                    .var$labels<-lf.contrastLabels(.var$levels,.var$contrast) 
                                    self$factors[[.var$name64]]<-.var
                                    self$names64$addFactor(.var$name,.var$levels)
                                    self$names64$addLabel(.var$name,.var$labels)
                                  }
                                  
                                  for (.var in self$covs) {
                                    self$names64$addVar(.var$name)
                                  }                          
                                },
                                 getInitData=function() {
                                   return(self$data64)
                                 },
                                cleandata=function(data) {
                                  data64<-list()
                                  data64[[self$dep$name64]]<-data[[self$dep$name]]
                                  
                                  for (i in seq_along(self$covs)) {
                                    .var<-self$covs[[i]]
                                    .var$oparams<-list(mean=mean(data[[.var$name]],na.rm = T),sd=sd(data[[.var$name]],na.rm = T))
                                    data64[[.var$name64]]<-private$.scaleContinuous(data[[.var$name]],method = .var$scale,by=self$clusters[[1]])
                                    attr(data64[[.var$name64]],"infos")<-.var
                                    self$covs[[i]]<-.var
                                  }
                                  for (.var in self$factors) {
                                    if (!is.factor(data[[.var$name]])) data[[.var$name]]<-factor(data[[.var$name]])
                                    data64[[.var$name64]]<-data[[.var$name]]
                                    stats::contrasts(self$data64[[.var$name64]]) <- lf.createContrasts(.var$levels, .var$contrast)
                                    colnames(stats::contrasts(data64[[.var$name64]]))<-.var$cnames64
                                    attr(data64[[.var$name64]],"infos")<-.var
                                  }
                                  attr(data64, 'row.names') <- seq_len(length(data64[[1]]))
                                  attr(data64, 'class') <- 'data.frame'      
                                  data64 <- jmvcore::naOmit(data64)
                                },
                                standardizedata=function(data) {
                                  zdata<-data
                                  zdata[[self$dep$name64]]<-as.numeric(scale(data[[self$dep$name64]]))
                                  for (.var in self$covs) {
                                    zdata[[.var$name64]]<-as.numeric(scale(data[[.var$name64]]))
                                  }
                                  zdata
                                },
                                
                                 print=function() {
                                   cat("###### covariates ######\n")
                                   lapply(self$covs,print)
                                   cat("###### factors ######\n")
                                    lapply(self$factors,print)
                                    cat("###### options ######\n")
                                    lapply(self$options,print)
                                    
                                 },
                                 formula=function() {
                                   dep<-self$dep$name
                                   fi<-self$options$fixedIntercept
                                   modelTerms<-self$options$modelTerms
                                   lf.constructFormula(dep,modelTerms,fi)
                                 },
                                 formula64=function() {
                                   dep64<-self$dep$name64
                                   fi<-self$options$fixedIntercept
                                   modelTerms64<-lapply(self$options$modelTerms,jmvcore::toB64)
                                   as.formula(lf.constructFormula(dep64,modelTerms64,fi))
                                },
                                 anovaTerms=function() {
                                   mynames64<-attr(terms(self$formula64()),"term.labels")
                                   self$names64$nicenames(mynames64)
                                 },
                                paramTerms=function() {
                                  mynames64<-colnames(model.matrix(self$formula64(),self$data64))
                                  self$names64$nicenames(mynames64)
                                },
                                
                                paramLabels=function() {
                                  mynames64<-colnames(model.matrix(self$formula64(),self$data64))
                                  self$names64$nicelabels(mynames64)
                                } 

                                ),
                            private=list(
                                .initOptions=function(options) {
                                  
                                  .options<-list()
                                  .options[["dep"]]<-options$dep
                                  .options[["covs"]]<-options$covs
                                  .options[["factors"]]<-options$factors

                                  if ("cluster" %in% names(options)) .options[["clusters"]]<-options$cluster else .options[["clusters"]]<-NULL                                  
                                  .contrasts<-lapply(options$contrasts,function(a) a$type)
                                  names(.contrasts)<-sapply(options$contrasts,function(a) a$var)
                                  .options[["contrasts"]]<-.contrasts
                                  .scaling<-sapply(options$scaling,function(a) a$type)
                                  names(.scaling)<-sapply(options$scaling,function(a) a$var)
                                  .options[["scaling"]]<-.scaling
                                  
                                  modelTerms<-options$modelTerms
                                  fixedIntercept<-options$fixedIntercept
                                  
                                  # this allows intercept only model to be passed by syntax interface
                                  aOne<-which(unlist(modelTerms)=="1")
                                  if (is.something(aOne)) {
                                    modelTerms[[aOne]]<-NULL
                                    fixedIntercept=TRUE
                                  }
                                  .options[["fixedIntercept"]]<-fixedIntercept
                                  .options[["modelTerms"]]<-modelTerms
                                   
                                  .type<-options$simpleScale
                                  .span<-ifelse(.type=="mean_sd",options$cvalue,options$percvalue)
                                  .options[["conditioning"]]<-list(type=.type,span=.span)
                                   self$options<-.options
                                } ,
                                .init=function() {
                                     options<-self$options
                                  ### here we gather all info about the variables
                                  
                                     ## dependent variable 
                                     
                                     self$dep<-list(name=self$options$dep,
                                                    name64=jmvcore::toB64(self$options$dep))

                                     #### continuous variables ######
                                     self$covs<-lapply(options$covs,function(a) {
                                         if (a %in% names(options$scaling)) .scale<-options$scaling[[a]] else .scale="centered"
                                         list(name=a,
                                              name64=jmvcore::toB64(a),
                                              scale=.scale)
                                      })
                                     names(self$covs)<-jmvcore::toB64(options$covs)
                                     
                                     ##### factors #######
                                     self$factors<-lapply(options$factors,function(a) {
                                       name<-a
                                       name64<-jmvcore::toB64(name)
                                       if (a %in% names(options$contrasts)) .type<-options$contrasts[[a]] else .type<-"simple"
                                        list(name=name,
                                             name64=jmvcore::toB64(name),
                                             contrast=.type
                                             )
                                     })
                                     names(self$factors)<-jmvcore::toB64(options$factors)
    
                                     if (is.something(options$clusters))
                                          self$clusters<-options$cluster
                                     self$vars<-unlist(c(options$dep,options$factors,options$cov,self$cluster))

                                },
                                .scaleContinuous=function(var,method,by=NULL) {
                                    var<-jmvcore::toNumeric(var)
                                  if (method=="centered") 
                                    var<-scale(var,scale = F)  
                                  if (method=="clusterbasedcentered") 
                                    var<-unlist(tapply(var,by,scale,scale=F))
                                  if (method=="standardized") 
                                    var<-scale(var,scale = T)  
                                  if (method=="clusterbasedstandardized")     
                                    var<-unlist(tapply(var,by,scale,scale=T))
                                  if (method=="log") {
                                    if (any(var<=0))
                                      jmvcore::reject("Log transformation requires all positive numbers in the original variable")
                                    else
                                      var<-log(var)
                                  }
                                  as.numeric(var)
                                }
                                
                                )
                                
                              
)

