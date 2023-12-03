# Репозитарий с кодом для IC SSA #

Необходимый список пакетов для работы:

1. [R](https://www.r-project.org/)
2. [Rtools](https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html) если используется Windows
3. Пакеты, необходимые для работы в R:
```
install.packages(c("svd", "Matrix", "Rssa", "devtools", "snow", "parallel"))
devtools::install_github("asl/ssabook")
devtools::install_github("furiousbean/rhlra", ref="ic_work_report")
```
