#' Merge a list or matrix of p-values
#'
#' @param scores Either a list of p-values or a matrix where each column is a test.
#' @param method Method to merge p-values. See 'methods' section below.
#'
#' @return If \code{scores} is a vector or list, returns a number. If \code{scores} is a
#'   matrix, returns a named list of p-values merged by row.
#'
#' @section Methods:
#' Two methods are available to merge a list of p-values:
#' \describe{
#'  \item{Fisher}{Fisher's method (default) assumes that p-values are uniformly
#'  distributed and performs a chi-squared test on the statistic sum(-2 log(p)).
#'  This method is most appropriate when the columns in \code{scores} are
#'  independent.}
#'  \item{Brown}{Brown's method extends Fisher's method by accounting for the
#'  covariance in the columns of \code{scores}. It is more appropriate when the
#'  tests of significance used to create the columns in \code{scores} are not
#'  necessarily independent. Note that the "Brown" method cannot be used with a 
#'  single list of p-values. However, in this case Brown's method is identical 
#'  to Fisher's method and should be used instead.}
#' }
#'
#' @examples
#'   merge_p_values(c(0.05, 0.09, 0.01))
#'   merge_p_values(list(a=0.01, b=1, c=0.0015, d=0.025), method='Fisher')
#'   merge_p_values(matrix(data=c(0.03, 0.061, 0.48, 0.052), nrow = 2), method='Brown')
#' 
#' @export
merge_p_values <- function(scores, method = "Fisher") {
    # Validation on scores
    if (is.list(scores)) scores <- unlist(scores, recursive=FALSE)
    if (!(is.vector(scores) || is.matrix(scores))) stop("scores must be a matrix or list.")
    if (any(is.na(scores))) stop("scores may not contain missing values.")
    if (!is.numeric(scores)) stop("scores must be numeric.")
    if (any(scores < 0 | scores > 1)) stop("All values in scores must be in [0,1].")
	if (!method %in% c("Fisher", "Brown")) stop("Only Fisher's and Brown's methods are currently supported.")

    if (is.vector(scores)) {

        if (method == "Brown") {
        	stop("Brown's method cannot be used with a single list of p-values")
        }
        
        # convert zeroes to smallest available doubles
        scores <- sapply(scores, function(x) ifelse (x == 0, 1e-300, x))
        
        # if not brown, then fisher
        return(fishersMethod(scores))
    }
    
    # scores is a matrix with one column, then no transformatino needed
    if (ncol(scores) == 1) return (scores[, 1, drop=TRUE])
    
    if (method == "Brown") {
        cov.matrix <- calculateCovariances(t(scores))
        return(apply(scores, 1, brownsMethod, cov.matrix=cov.matrix))
    }
    
    scores <- apply(scores, c(1,2), function(x) ifelse (x == 0, 1e-300, x))

    return (apply(scores, 1, fishersMethod))
}


fishersMethod <- function(p.values) {
    lnp <- log(p.values)
    chisq <- (-2) * sum(lnp)
    df <- 2 * length(lnp)
    stats::pchisq(chisq, df, lower.tail = FALSE)
}


#' Merge p-values using the Brown's method. 
#'
#' @param p.values A vector of m p-values.
#' @param data.matrix An m x n matrix representing m tests and n samples. NA's are not allowed.
#' @param cov.matrix A pre-calculated covariance matrix of \code{data.matrix}. This is more
#'   efficient when making many calls with the same data.matrix.
#'   Only one of \code{data.matrix} and \code{cov.matrix} must be given. If both are supplied,
#'   \code{data.matrix} is ignored.
#' @return A single p-value representing the merged significance of multiple p-values.
#' @export

# Based on the R package EmpiricalBrownsMethod
# https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM/blob/master/R/EmpiricalBrownsMethod/R/ebm.R
# Only significant differences are the removal of extra_info and allowing a
# pre-calculated covariance matrix
# 
brownsMethod <- function(p.values, data.matrix = NULL, cov.matrix = NULL) {
    if (missing(data.matrix) && missing(cov.matrix)) {
        stop ("Either data.matrix or cov.matrix must be supplied")
    }
    if (!(missing(data.matrix) || missing(cov.matrix))) {
        message("Both data.matrix and cov.matrix were supplied. Ignoring data.matrix")
    }
    if (missing(cov.matrix)) cov.matrix <- calculateCovariances(data.matrix)

    N <- ncol(cov.matrix)
    expected <- 2 * N
    cov.sum <- 2 * sum(cov.matrix[lower.tri(cov.matrix, diag=FALSE)])
    var <- (4 * N) + cov.sum
    sf <- var / (2 * expected)

    df <- (2 * expected^2) / var
    if (df > 2 * N) {
        df <- 2 * N
        sf <- 1
    }

    x <- 2 * sum(-log(p.values), na.rm=TRUE)
    p.brown <- stats::pchisq(df=df, q=x/sf, lower.tail=FALSE)
    p.brown
}

transformData <- function(dat) {
    # If all values in dat are the same (equal to y), return dat. The covariance
    # matrix will be the zero matrix, and brown's method gives the p-value as y
    # Otherwise (dat - dmv) / dvsd is NaN and ecdf throws and error
    if (isTRUE(all.equal(min(dat), max(dat)))) return(dat)

    dvm <- mean(dat, na.rm=TRUE)
    dvsd <- pop.sd(dat)
    s <- (dat - dvm) / dvsd
    distr <- stats::ecdf(s)
    sapply(s, function(a) -2 * log(distr(a)))
}


calculateCovariances <- function(data.matrix) {
    transformed.data.matrix <- apply(data.matrix, 1, transformData)
    stats::cov(transformed.data.matrix)
}

pop.var <- function(x) stats::var(x, na.rm=TRUE) * (length(x) - 1) / length(x)
pop.sd <- function(x) sqrt(pop.var(x))



#'
#' Aggregated Cauchy Assocaition Test
#'
#' A p-value combination method using the Cauchy distribution.
#'
#'
#'
#' @param weights a numeric vector/matrix of non-negative weights for the combined p-values. When it is NULL, the equal weights are used.
#' @param Pvals a numeric vector/matrix of p-values. When it is a matrix, each column of p-values is combined by ACAT.
#' @param is.check logical. Should the validity of \emph{Pvals} (and \emph{weights}) be checked? When the size of \emph{Pvals} is large and one knows \emph{Pvals} is valid, then the checking part can be skipped to save memory.
#' @return The p-value(s) of ACAT.
#' @author Yaowu Liu
#' @examples p.values<-c(2e-02,4e-04,0.2,0.1,0.8,04e-10,4e-10,4e-04);ACAT(Pvals=p.values)
#' @examples ACAT(matrix(runif(1000),ncol=10))
#' @references Liu, Y., & Xie, J. (2019). Cauchy combination test: a powerful test with analytic p-value calculation
#' under arbitrary dependency structures. \emph{Journal of American Statistical Association},115(529), 393-402. (\href{https://amstat.tandfonline.com/doi/abs/10.1080/01621459.2018.1554485}{pub})
#' @export
#' 

ACAT<-function(Pvals,weights=NULL,is.check=TRUE){
    Pvals<-as.matrix(Pvals)
    if (is.check){
        #### check if there is NA
        if (sum(is.na(Pvals))>0){
            stop("Cannot have NAs in the p-values!")
        }
        #### check if Pvals are between 0 and 1
        if ((sum(Pvals<0)+sum(Pvals>1))>0){
            stop("P-values must be between 0 and 1!")
        }
        #### check if there are pvals that are either exactly 0 or 1.
        is.zero<-(colSums(Pvals==0)>=1)
        is.one<-(colSums(Pvals==1)>=1)
        if (sum((is.zero+is.one)==2)>0){
            stop("Cannot have both 0 and 1 p-values in the same column!")
        }

        if (sum(is.zero)>0){
            warning("There are p-values that are exactly 0!")
        }
        if (sum(is.one)>0){
            warning("There are p-values that are exactly 1!")
        }

    }
    #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
    if (is.null(weights)){
        is.weights.null<-TRUE
    }else{
        is.weights.null<-FALSE
        weights<-as.matrix(weights)
        if (sum(dim(weights)!=dim(Pvals))>0){
            stop("The dimensions of weights and Pvals must be the same!")
        }else if (is.check & (sum(weights<0)>0)){
            stop("All the weights must be nonnegative!")
        }else{
            w.sum<-colSums(weights)
            if (sum(w.sum<=0)>0){
                stop("At least one weight should be positive in each column!")
            }else{
                for (j in 1:ncol(weights)){
                    weights[,j]<-weights[,j]/w.sum[j]
                }
            }
        }

    }

    #### check if there are very small non-zero p values and calcuate the cauchy statistics
    is.small<-(Pvals<1e-15)
    if (is.weights.null){
         Pvals[!is.small]<-tan((0.5-Pvals[!is.small])*pi)
         Pvals[is.small]<-1/Pvals[is.small]/pi
         cct.stat<-colMeans(Pvals)
    }else{
         Pvals[!is.small]<-weights[!is.small]*tan((0.5-Pvals[!is.small])*pi)
         Pvals[is.small]<-(weights[is.small]/Pvals[is.small])/pi
         cct.stat<-colSums(Pvals)
    }
    #### return the ACAT p value(s).
    pval<-pcauchy(cct.stat,lower.tail = F)
    return(pval)
}

#'
#' A set-based test that uses ACAT to combine the variant-level p-values.
#'
#'
#' @param G a numeric matrix or dgCMatrix with each row as a different individual and each column as a separate gene/snp. Each genotype should be coded as 0, 1, 2.
#' @param obj an output object of the \code{\link{NULL_Model}} function.
#' @param weights.beta a numeric vector of parameters for the beta weights for the weighted kernels. If you want to use your own weights, please use the “weights” parameter. It will be ignored if “weights” parameter is not null.
#' @param weights a numeric vector of weights for the SNP p-values. When it is NULL, the beta weight with the “weights.beta” parameter is used.
#' @param mac.thresh a threshold of the minor allele count (MAC). The Burden test will be used to aggregate the SNPs with MAC less than this thrshold.
#' @return The p-value of ACAT-V.
#' @details The Burden test is first used to aggregate very rare variants with Minor Allele Count (MAC) < \emph{mac.thresh} (e.g., 10), and a Burden p-value is obtained. For each of the variants with MAC >= \emph{mac.thresh}, a variant-level p-value is calculated. Then, ACAT is used to combine the variant-level p-values and the Burden test p-value of very rare variants.
#'
#' If \emph{weights.beta} is used, then the weight for the Burden test p-value is demetermined by the average Minor Allele Frequency (MAF) of the variants with MAC < \emph{mac.thresh}; if the user-specified \emph{weights} is used, then the weight for the Burden test p-value is the average of \emph{weights} of the variants with MAC < \emph{mac.thresh}.
#'
#' Note that the \emph{weights} here are for the SNP p-vlaues. In SKAT, the weights are for the SNP score test statistics. To transfrom the SKAT weights to the \emph{weights} here, one can use the formula that \emph{weights} = (skat_weights)^2*MAF*(1-MAF).
#'
#' @author Yaowu Liu
#' @references Liu, Y., et al. (2019). ACAT: A fast and powerful p value combination
#' method for rare-variant analysis in sequencing studies.
#' \emph{American Journal of Humann Genetics 104}(3), 410-421.
#' (\href{https://www.sciencedirect.com/science/article/pii/S0002929719300023}{pub})
#' @examples  library(Matrix)
#' @examples  data(Geno)
#' @examples  G<-Geno[,1:100] # Geno is a dgCMatrix of genotypes
#' @examples  Y<-rnorm(nrow(G)); Z<-matrix(rnorm(nrow(G)*4),ncol=4)
#' @examples  obj<-NULL_Model(Y,Z)
#' @examples  ACAT_V(G,obj)
#' @export
ACAT_V<-function(G,obj,weights.beta=c(1,25),weights=NULL,mac.thresh=10){
    ### check weights
    if (!is.null(weights) && length(weights)!=ncol(G)){
        stop("The length of weights must equal to the number of variants!")
    }

    mac<-Matrix::colSums(G)
    ### remove SNPs with mac=0
    if (sum(mac==0)>0){
        G<-G[,mac>0,drop=FALSE]
        weights<-weights[mac>0]
        mac<-mac[mac>0]
        if (length(mac)==0){
            stop("The genotype matrix do not have non-zero element!")
        }
    }
    ### p and n
    p<-length(mac)
    n<-nrow(G)
    ###

    if (sum(mac>mac.thresh)==0){  ## only Burden
        pval<-Burden(G,obj, weights.beta = weights.beta, weights = weights)
    }else if (sum(mac<=mac.thresh)==0){ ## only cauchy method
        if (is.null(weights)){
            MAF<-mac/(2*n)
            W<-(dbeta(MAF,weights.beta[1],weights.beta[2])/dbeta(MAF,0.5,0.5))^2
        }else{
            W<-weights
        }

        Mpvals<-Get.marginal.pval(G,obj)
        pval<-ACAT(Mpvals,W)
    }else{  ## Burden + Cauchy method
        is.very.rare<-mac<=mac.thresh
        weights.sparse<-weights[is.very.rare]
        weights.dense<-weights[!is.very.rare]
        pval.dense<-Burden(G[,is.very.rare,drop=FALSE],obj, weights.beta = weights.beta, weights = weights.sparse)

        Mpvals<-Get.marginal.pval(G[,!is.very.rare,drop=FALSE],obj)

        Mpvals<-c(Mpvals,pval.dense)
        if (is.null(weights)){
            MAF<-mac/(2*n)
            mafs<-c(MAF[!is.very.rare],mean(MAF[is.very.rare])) ## maf for p-values
            W<-(dbeta(mafs,weights.beta[1],weights.beta[2])/dbeta(mafs,0.5,0.5))^2
        }else{
            W<-c(weights.dense,mean(weights.sparse))
        }


        is.keep<-rep(T,length(Mpvals))
        is.keep[which(Mpvals==1)]<-F  ## remove p-values of 1.
        pval<-ACAT(Mpvals[is.keep],W[is.keep])
    }
    return(pval)
}

#'
#'
#' Get parameters and residuals from the NULL model
#'
#' Compute model parameters and residuals for ACAT-V
#'
#'
#' @param Y a numeric vector of outcome phenotypes.
#' @param Z a numeric matrix of covariates. Z must be full-rank. Do not include intercept in Z. The intercept will be added automatically.
#' @return This function returns an object that has model parameters and residuals of the NULL model of no association between genetic variables and outcome phenotypes. After obtaining it, please use \code{\link{ACAT_V}} function to conduct the association test.
#' @details \emph{Y} could only be continuous or binary. If \emph{Y} is continuous, a linear regression model is fitted. If \emph{Y} is binary, it must be coded as 0,1 and a logistic model is fitted.
#' @author Yaowu Liu
#' @examples  Y<-rnorm(10000)
#' @examples  Z<-matrix(rnorm(10000*4),ncol=4)
#' @examples  obj<-NULL_Model(Y,Z)
#' @export
NULL_Model<-function(Y,Z=NULL){
    n<-length(Y)
    #### check the type of Y
    if ((sum(Y==0)+sum(Y==1))==n){
        out_type<-"D"
    }else{
        out_type<-"C"
    }
    #### Add intercept
    Z.tilde<-cbind(rep(1,length(Y)),Z)
    if (out_type=="C"){
        #### estimate of sigma square
        Z.med<-Z.tilde%*%solve(chol(t(Z.tilde)%*%Z.tilde))   ## Z.med%*%t(Z.med) is the projection matrix of Z.tilde
        Y.res<-as.vector(Y-(Y%*%Z.med)%*%t(Z.med))
        sigma2<-sum(Y.res^2)/(n-ncol(Z.med))
        #### output
        res<-list()
        res[["out_type"]]<-out_type
        res[["Z.med"]]<-Z.med
        res[["Y.res"]]<-Y.res
        res[["sigma2"]]<-sigma2
    }else if (out_type=="D"){
        #### fit null model
        g<-glm(Y~0+Z.tilde,family = "binomial")
        prob.est<-g[["fitted.values"]]
        #### unstandarized residuals
        Y.res<-(Y-prob.est)
        ### Sigma when rho=0
        sigma2.Y<-prob.est*(1-prob.est)  ### variance of each Y_i
        ### output
        res<-list()
        res[["out_type"]]<-out_type
        res[["Z.tilde"]]<-Z.tilde
        res[["Y.res"]]<-Y.res
        res[["sigma2.Y"]]<-sigma2.Y
    }
    return(res)
}




Get.marginal.pval<-function(G,obj){
    ### check obj
    if (names(obj)[1]!="out_type"){
        stop("obj is not calculated from MOAT_NULL_MODEL!")
    }else{
        out_type<-obj[["out_type"]]
        if (out_type=="C"){
            if (!all.equal(names(obj)[2:length(obj)],c("Z.med","Y.res","sigma2"))){
                stop("obj is not calculated from MOAT_NULL_MODEL!")
            }else{
                Z.med<-obj[["Z.med"]]
                Y.res<-obj[["Y.res"]]
                n<-length(Y.res)
                SST<-obj[["sigma2"]]*(n-ncol(Z.med))
            }
        }else if (out_type=="D"){
            if (!all.equal(names(obj)[2:length(obj)],c("Z.tilde","Y.res","sigma2.Y"))){
                stop("obj is not calculated from MOAT_NULL_MODEL!")
            }else{
                Z.tilde<-obj[["Z.tilde"]]
                Y.res<-obj[["Y.res"]]
                sigma2.Y<-obj[["sigma2.Y"]]
                n<-length(Y.res)
            }
        }
    }

    if (class(G)!="matrix" && class(G)!="dgCMatrix"){
        stop("The class of G must be matrix or dgCMatrix!")
    }

    if (out_type=="C"){
        G_tX.med<-as.matrix(Matrix::crossprod(Z.med,G))
        ### Sigma^2 of G
        Sigma2.G<-Matrix::colSums(G^2)-Matrix::colSums(G_tX.med^2)
        SSR<-as.vector((Y.res%*%G)^2/Sigma2.G)
        SSR[Sigma2.G<=0]<-0
        df.2<-n-1-ncol(Z.med)
        t.stat<-suppressWarnings(sqrt(SSR/((SST-SSR)/df.2)))
        marginal.pvals<-2*pt(t.stat,(n-1-ncol(Z.med)),lower.tail = FALSE)
    }else if (out_type=="D"){
        Z.stat0<-as.vector(Y.res%*%G)
        ### Sigma when rho=0
        tG_X.tilde_sigma2<-as.matrix(Matrix::crossprod(G,Z.tilde*sigma2.Y))
        Sigma2.G<-Matrix::colSums(G^2*sigma2.Y)-diag(tG_X.tilde_sigma2%*%solve(t(Z.tilde)%*%(Z.tilde*sigma2.Y))%*%t(tG_X.tilde_sigma2))
        marginal.pvals<-2*pnorm(abs(Z.stat0)/sqrt(Sigma2.G),lower.tail = FALSE)
    }

    return(marginal.pvals)
}


Burden<-function(G,obj,kernel="linear.weighted",weights.beta=c(1,25),weights=NULL){
    ### check obj
    if (names(obj)[1]!="out_type"){
        stop("obj is not calculated from NULL_MODEL!")
    }else{
        out_type<-obj[["out_type"]]
        if (out_type=="C"){
            if (!all.equal(names(obj)[2:length(obj)],c("Z.med","Y.res","sigma2"))){
                stop("obj is not calculated from NULL_MODEL!")
            }else{
                Z.med<-obj[["Z.med"]]
                Y.res<-obj[["Y.res"]]/sqrt(obj[["sigma2"]])  ## rescaled residules
                n<-length(Y.res)
            }
        }else if (out_type=="D"){
            if (!all.equal(names(obj)[2:length(obj)],c("Z.tilde","Y.res","sigma2.Y"))){
                stop("obj is not calculated from NULL_MODEL!")
            }else{
                Z.tilde<-obj[["Z.tilde"]]
                Y.res<-obj[["Y.res"]]
                sigma2.Y<-obj[["sigma2.Y"]]
                n<-length(Y.res)
            }
        }
    }
    ### MAF
    MAF<-Matrix::colSums(G)/(2*dim(G)[1])
    p<-length(MAF)
    #### weights
    if (kernel=="linear.weighted"){
        if (is.null(weights)){
            W<-dbeta(MAF,weights.beta[1],weights.beta[2])
        }else{
            if (length(weights)==p){
                W<-weights
            }else{
                stop("The length of weights must equal to the number of variants!")
            }
        }

    }else if (kernel=="linear"){
        W<-rep(1,p)
    }else{
        stop("The kernel name is not valid!")
    }

    ###### if G is sparse or not
    if (class(G)=="matrix" || class(G)=="dgCMatrix"){
        if (out_type=="C"){
            Z.stat.sum<-as.vector((Y.res%*%G)%*%W)
            Gw<-G%*%W
            sigma.z<-sqrt(sum(Gw^2)-sum((t(Z.med)%*%(Gw))^2))
        }else if (out_type=="D"){
            Z.stat.sum<-as.vector((Y.res%*%G)%*%W)
            Gw<-as.vector(G%*%W)
            sigma.z<-sum(Gw^2*sigma2.Y)-((Gw*sigma2.Y)%*%Z.tilde)%*%solve(t(Z.tilde)%*%(Z.tilde*sigma2.Y))%*%t((Gw*sigma2.Y)%*%Z.tilde)
            sigma.z<-as.vector(sqrt(sigma.z))
        }
    }else{
        stop("The class of G must be matrix or dgCMatrix!")
    }

    V<-Z.stat.sum/sigma.z
    Q<-V^2   ## Q test statistic
    pval<-1-pchisq(Q,df=1)
    return(pval)
}