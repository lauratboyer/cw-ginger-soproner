
# Opposite of %in%
"%nin%" <- function(x, y) !x %in% y
"%val.in%" <- function(x, y) x[x %in% y]
"%val.nin%" <- function(x, y) x[!x %in% y]

# Get the number of unique values
count <- function(x) length(unique(x))

# Timer functions
start.timer <- function() assign("timer",proc.time()[3],.GlobalEnv)
stop.timer <- function() print(proc.time()[3]-timer)

# Looks for a pattern in objects named in the Global environment
objFind <- function(x) ls(.GlobalEnv)[grep(x,ls(.GlobalEnv),ignore.case=TRUE)]

# Handle for grep(..., value=TRUE)
grepv <- function(...) grep(..., value=TRUE)

####################################################
# table2pdf()
# Creates pdf table from matrix or dataframe

# Note 1: # You need a working version of latex in your computer. Make sure
# latex in your R search path. Sys.getenv("PATH") will show you the value
# of your current path, see ?Sys.getenv and Details section in ?latex for more.
if(!grepl("tex",Sys.getenv("PATH"))) {
    print("slow down cowboy! latex is either
           not installed or not on your path")}

# Note 2: "capt" is by default set at NULL, which results in no heading to the table.
# To have a simple header with the table number but no caption, set capt=" "
# The "num" argument controls the table number in the table caption, with
# num=0 producing a simple "Table: " header with no numbers

# Note 3: crops to remove white margins, could be set up
# otherwise if pdfcrop function is not available (Hello there, PCs!)
# could use latex package anysize to set paper size and margins
# as a function of number of rows/columns and number of characters in each
# Apparently Ghostscript/Ghostview will do the cropping for PCs
# see: http://www.ats.ucla.edu/stat/latex/icu/install_win.htm
table2pdf <- function(data, capt=NULL, num=1, file="temptable.pdf",
                      useRowNames=FALSE,rownames.col,keepTex=FALSE,...) {

# Sample table if no data given in argument:
if(missing(data))
  {data <- matrix(1:6, nrow=2, dimnames=list(c('a','b'),c('c','d','this that')))
   capt <- "A sample table just for fun"}

# Creating formatted .tex file
preamble <- list(
                 paste("\\","documentclass{report} \n",sep=""),
                 paste("\\","usepackage{booktabs} \n",sep=""),
                 paste("\\","usepackage{pdflscape} \n",sep=""),
                 paste("\\","usepackage{color,colortbl} \n",sep=""),
                 paste("\\","thispagestyle{empty} \n",sep=""),
                 paste("\\","begin{document} \n",sep=""),
                 paste("\\","setcounter{table}{",num-1,"} \n", sep=""),
                 paste("\\","definecolor{Gray}{gray}{0.9} \n", sep=""),
                 paste("\\","definecolor{LightCyan}{rgb}{0.88,1,1}", sep="")

                 )

if(num==0) preamble[[length(preamble)+1]] <- paste("\\","renewcommand","\\","thetable{} \n",sep="")

cat(paste(preamble),file="temptable.tex")

rl <- "data"
if(!missing(rownames.col)) rl <- rownames.col
if(!useRowNames) {
    rn <- NULL
    rl <- NULL
    #rownames(data) <- data[,1]
    #              rl <- names(data)[1]
              }else{
                  rn <- rownames(data)}

# Format column names so that latex likes them
names(data) <- gsub("_","\\\\_",names(data))

latex(data,file="temptable.tex",append=TRUE, caption=capt, booktabs=TRUE,
      rowlabel=rl,rowname=rn,...)
cat(paste("\\","end{document}",sep=""),file="temptable.tex", append=TRUE)

# Convert .tex to .pdf
texi2dvi("temptable.tex",pdf=T)
system(paste("pdfcrop temptable.pdf ", "temptable.pdf", " --margins 5", sep=""),
       intern=TRUE)
file.rename("temptable.pdf",file)
if(!keepTex) file.remove("temptable.tex")
fr <- file.remove(c("temptable.aux","temptable.log"))
}


## -------------------------------------------------------
## file-scan.r
## Function to find a pattern in files with user-defined extension (.r by default)
## Returns file name and approximate row matches if positive.


file.scan <- function(pattern, dir=getwd(), fpattern=NA, ftype=".r$|.R$") {

    nfiles <- list.files(dir,full.names=TRUE)
    if(!is.na(fpattern)) nfiles <- nfiles[grep(fpattern, nfiles)]
    if(!(ftype=="all")) nfiles <- nfiles[grep(sprintf("%s$",ftype),nfiles)]
    if(length(nfiles)==0) print("No files found")


    sfunk <- function(wf) {
        ss <- scan(wf, "character", sep="\n", blank.lines.skip=FALSE,quiet=TRUE)
        gf <- grep(pattern, ss)
        if(length(gf)>0) {
            return(list("filename"=wf, "row.numbers"=gf))
        }else{return(NULL)}
    }

    fs <- sapply(nfiles, sfunk)
    fs <- fs[!sapply(fs,is.null)]

    return(names(fs))
}
