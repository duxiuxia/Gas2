detach(package:JanusPackge)
remove.packages(pkgs="JanusPackge")
rm(list=ls())
package.skeleton(name="JanusPackge",path=".",force=TRUE,namespace=TRUE,code_files=c("code/pipeline.r"))
R CMD check Janus/JanusPackge
R CMD build Janus/JanusPackge --binary
R CMD install JanusPackge_1.0_R_i386-apple-darwin8.11.1.tar.gz

Rscript Janus/code/webUI.r /Users/mike/Janus/workDir
zip -m "TIC.zip" ./*
