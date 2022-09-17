

#---------------------------- Function for diagnostic analysis -------------------------------#
DiagFunc <-
  function(Rep = "REP",
           Trait = NULL,
           Trat = NULL,
           data1 = data,
           Exp = "Test",
           plot_diag1 = TRUE,
           plot_diag2 = TRUE,
           plotBox1 = TRUE,
           lambda.boxcox = c(0.9, 1.1),
           # Número de análises realizadas simultaneamente
           nDiag = 3,
           verbose = FALSE,
           # Column names to return
           ColumnNames_To_Return =
             c(NULL, NULL)) {
    # Input Exp column if necessary
    if (is.null(Exp)) {
      data1$Exp <- "Exp1"
      Exp <- "Exp"
      Exp2 <- NULL
    } else {
      Exp2 <- "Exp2"
    }
    
    
    if (!require("pacman")) {
      install.packages("pacman")
    }
    
    pacman::p_load(MASS, lmerTest, nortest, car, moments, Hmisc, tidyverse)
    
    ColumnNames_To_Return <- c(ColumnNames_To_Return, Trait)
    
    if (is.null(Rep)) {
      data1$REP <- 1
      Rep <- "REP"
    }
      colnames(data1)[match(c(Rep, Trat, Exp, Trait), colnames(data1))] <-
        c("Rep", "Trat", "Exp", "Trait")
    
    
    if (length(colnames(data1)[duplicated(colnames(data1))]) != 0) {
      stop(paste("Please, check and rename the following columns in the dataset:",
                 colnames(data1)[duplicated(colnames(data1))]))
    }
    
    # Test and boxcox transformation, if necessary
    FIT <- list()
    BOX <- list()
    LAMBDA <- list()
    DATA <- list()
    MODEL <- list()
    for (i in unique(unlist(data1[["Exp"]]))) {
      if (!is.null(Rep) & (
        data1 %>%
        dplyr::filter(Exp == "exp1") |>
        dplyr::select(Rep, Trait) |>
        filter(!is.na(Trait) & !duplicated(Rep)) |>
        dplyr::select(Trait) |>
        unlist() %>%
        sum(.) %>%
        is.na(.)
      ) == T) {
        Model <-
          as.formula(Trait ~ as.factor(Trat) + as.factor(Rep),
                     env = parent.frame(0))
      } else {
        Model <- as.formula(Trait ~ as.factor(Trat), env = parent.frame(0))
      }
      
      Model2 <-
        as.formula(Trait ~ as.factor(Trat), env = parent.frame(0))
      fit <-
        lm(Model,
           data = subset(data1, Exp == i))
      
      FIT[[i]] <- fit
      MODEL[[i]] <- Model
      DATA[[i]] <- subset(data1, Exp == i)
    }
    for (i in unique(unlist(data1[["Exp"]]))) {
      # boxcox
      box <-
        boxcox(
          MODEL[[i]],
          lambda = seq(-3, 3, 0.05),
          plotit = F,
          data = DATA[[i]]
        )
      lambda <- box$x[which(box$y == max(box$y))]
      
      BOX[[i]] <- box
      LAMBDA[[i]] <- lambda
      
    }
    
    
    for (i in unique(unlist(data1[["Exp"]]))) {
      DATA[[i]]$TraitT <- NA
      if (BOX[[i]]$x[which(BOX[[i]]$y == max(BOX[[i]]$y))] < lambda.boxcox[1] |
          BOX[[i]]$x[which(BOX[[i]]$y == max(BOX[[i]]$y))] > lambda.boxcox[2]) {
        cat(paste0(
          "\nTransformation to ",
          i,
          ": ",
          "lambda = ",
          round(LAMBDA[[i]], 3),
          "\n"
        ))
        # Transformando caractere (Trait)
        
        DATA[[i]]$TraitT <- DATA[[i]][["Trait"]] ^ LAMBDA[[i]]
        # data1 <- bind_rows(DATA)
        
        T_BoxCox <- TRUE
      } else {
        cat(paste0("\nNo transformation is required for ", i, "\n"))
        DATA[[i]]$TraitT <- DATA[[i]][["Trait"]]
        #data1 <- bind_rows(DATA)
        T_BoxCox <- FALSE
      }
    }
    data1 <- bind_rows(DATA)
    rm(DATA)
    
    # if (plotBox1 == TRUE) {
    #   for (i in unique(unlist(data1[["Exp"]]))) {
    #     boxcox(
    #       FIT[[i]],
    #       seq(-3, 3, 0.1),
    #
    #       xlab = bquote( ~ lambda == .(LAMBDA[[i]]))
    #     )
    #     title(i)
    #   }
    # }
    data1 <- data1
    
    colnames(data1)[match(c("Rep", "Trat", "Exp", "Trait"), colnames(data1))] <-
      c(Rep, Trat, Exp, Trait)
    
    if (is.null(Exp2)) {
      data1 <- data1 |>
        rename(Experiment = Exp)
      
      Exp <- "Experiment"
    }
    
    if (!is.null(Rep) & (
      data1 %>%
      filter(get(Exp) == i) |>
      dplyr::select(Rep, Trait) |>
      filter(!is.na(Trait) & !duplicated(get(Rep))) |>
      dplyr::select(Trait) |>
      unlist() %>%
      sum(.) %>%
      is.na(.)
    ) == F & unique(data1[[Rep]]) > 1) {
      Model2 <-
        as.formula(TraitT ~ as.factor(get(Trat)) + as.factor(get(Rep)))
    } else {
      Model2 <- as.formula(TraitT ~ as.factor(get(Trat)))
    }
    
    MOD <- list()
    for (i in unique(unlist(data1[[Exp]]))) {
      mod1 <-
        lm(Model2,
           data = data1 %>%
             filter(get(Exp) == i),
           na.action = na.omit)
      
      # mod2 = sem NA nos resíduos
      mod2 <-
        lm(Model2,
           data = data1 %>%
             filter(get(Exp) == i),
           na.action = na.exclude)
      MOD[[i]] <- list(
        mod1 = mod1,
        mod2 = mod2,
        DATA = data1 %>%
          filter(get(Exp) == i)
      )
    }
    
    
    for (i in unique(unlist(data1[[Exp]]))) {
      # Extraindo o resíduo estudentizado
      MOD[[i]]$DATA[, paste0("residuals.", Trait)] <-
        rstudent(MOD[[i]]$mod2, type = "pearson")
      # Extraindo os valores preditos
      MOD[[i]]$DATA[, paste0("predict.", Trait)] <-
        predict(MOD[[i]]$mod2)
      # Obtendo os resíduos estudentizados
      # Sem NA's
      MOD[[i]]$rs <- rstudent(MOD[[i]]$mod1, type = "pearson")
      
      # Com NA's
      MOD[[i]]$rs.n.na <- rstudent(MOD[[i]]$mod2, type = "pearson")
    }
    
    if (plot_diag1 == TRUE) {
      YP <- list()
      LMAX_RS <- list()
      LMIN_RS <- list()
      LMAX_YP <- list()
      LMIN_YP <- list()
      for (i in unique(unlist(data1[[Exp]]))) {
        if (max(MOD[[i]]$rs, na.rm = T) > 3 |
            min(MOD[[i]]$rs, na.rm = T) < -3) {
          DI1 <- "\nBefore removal of discrepant data"
        } else {
          DI1 <- ""
        }
      }
      # Escrevendo subtítulo nos gráficos
      SUBt <- list()
      for (i in unique(unlist(data1[[Exp]]))) {
        if (T_BoxCox == T) {
          SUBt[[i]] <-
            bquote(`BoxCox Transformation:` ~ lambda == .(LAMBDA[[i]]))
        } else {
          SUBt[[i]] <- ""
        }
      }
      for (i in unique(unlist(data1[[Exp]]))) {
        # Histograma
        hist(
          MOD[[i]]$rs.n.na,
          main = "",
          xlab = "Studentized Residuals",
          ylab = paste("Probability density"),
          freq = F
        )
        title(
          main = paste0("Histogram for ", Trait, " trait:", i, DI1),
          cex.main = 1,
          sub = SUBt[[i]]
        )
        xm <-
          seq(min(MOD[[i]]$rs, na.rm = T), max(MOD[[i]]$rs, na.rm = T), length = 40)
        ym <- dnorm(xm)
        lines(xm, ym)
        
        # qqNorm
        ## Plotando os pontos
        qqnorm(MOD[[i]]$rs.n.na,
               main = "")
        # Plotando a reta
        qqline(MOD[[i]]$rs.n.na,
               main = "",
               xlab = "",
               ylab = "")
        title(
          main = paste("Normality Graph - qqNorm:", Trait, "-", i, DI1),
          cex.main = 1,
          sub = SUBt[[i]]
        )
        # qqPlot
        qqPlot(MOD[[i]]$rs.n.na,
               main = "",
               ylab = "Studentized Residuals")
        title(
          main = paste("Normality Graph - qqPlot:", Trait, "-", i, DI1),
          cex.main = 1,
          sub = SUBt[[i]]
        )
        
        # Checando presença de valores discrepantes
        # Armazenando yp = Valores preditos
        YP[[i]] <- predict(MOD[[i]]$mod1, na.rm = T)
        
        # Marcando limites para valores normais
        LMAX_RS[[i]] <- max(MOD[[i]]$rs, 3.0, na.rm = T) + 0.1
        LMIN_RS[[i]] <- min(MOD[[i]]$rs, -3.0, na.rm = T) - 0.1
        
        # Armazenando os valores preditos
        LMAX_YP[[i]] <- max(YP[[i]], na.rm = T) + 0.1
        LMIN_YP[[i]] <- min(YP[[i]], na.rm = T) - 0.1
        
        # Plotando o gráfico
        ## Configurando os eixos e títulos
        plot(
          c(LMIN_YP[[i]], LMAX_YP[[i]]),
          c(LMIN_RS[[i]], LMAX_RS[[i]]),
          "n",
          main = "",
          xlab = expression(paste(hat(y)[predito])),
          ylab = "Studentized Residuals"
        )
        title(
          main = paste("Residual Analysis:", Trait, "-", i, DI1),
          cex.main = 1,
          sub = SUBt[[i]]
        )
        
        ## Plotando as linhas
        abline(h = 0, col = "black")
        abline(h = 3.0, col = "red")
        abline(h = -3.0, col = "red")
        
        ## Plotando os pontos
        points(YP[[i]], MOD[[i]]$rs)
      }
    } else {
      DI1 <- ""
      SUBt <- list()
      for (i in unique(unlist(data1[[Exp]]))) {
      SUBt[[i]] <- ""
      message("plot_diag1 = FALSE: omitting graphs before discrepant analysis")
      }
    }
    
    NORMAL1 <- list()
    ASSIMETRIA1 <- list()
    CURTOSE1 <- list()
    HOMOCED1 <- list()
    # Normalidade e homocedasticidade
    for (i in unique(unlist(data1[[Exp]]))) {
      # Teste de normalidade
      Normal <- list(ad.test(MOD[[i]]$rs), lillie.test(MOD[[i]]$rs))
      NORMAL1[[i]] <- Normal
      # Assimetria
      Assimetria <- skewness(MOD[[i]]$rs.n.na, na.rm = T)
      ASSIMETRIA1[[i]] <- Assimetria
      # Curtose
      Curtose <- kurtosis(MOD[[i]]$rs.n.na, na.rm = T)
      CURTOSE1[[i]] <- Curtose
      
      # Homocedasticidade
      homoced <- leveneTest(TraitT ~ as.factor(get(Trat)),
                            data = data1)
      HOMOCED1[[i]] <- homoced
      
      if (T_BoxCox == TRUE) {
        label(NORMAL1[[i]]) <-
          paste("\nNormality Test for",
                i,
                "after BoxCox transformation",
                DI1)
        label(ASSIMETRIA1[[i]]) <-
          paste("\nSkewness in the data of ",
                i,
                " after BoxCox transformation",
                DI1)
        label(CURTOSE1[[i]]) <-
          paste("\nKurtosis in the data of ",
                i,
                " after BoxCox transformation",
                DI1)
        label(HOMOCED1[[i]]) <-
          paste("\nHomoscedasticity Test for",
                i,
                " after BoxCox transformation",
                DI1)
      } else {
        label(NORMAL1[[i]]) <- paste("\nHomoscedasticity Test for", i, DI1)
        label(ASSIMETRIA1[[i]]) <-
          paste("\nSkewness in the data of", i, DI1)
        label(CURTOSE1[[i]]) <-
          paste("\nKurtosis in the data of", i, DI1)
        label(HOMOCED1[[i]]) <-
          paste("\nHomoscedasticity Test for", i, DI1)
      }
    }
    
    # Retirada de dados Discrepantes
    
    disc <- 0
    # rm(disc)
    # disc3 <- 0
    data2 <- data1
    # MOD2 <- list()
    
    repeat {
      data2 <- data2
      # MOD2 <- MOD
      
      if (!is.null(Rep) & (
        data1 %>%
        filter(get(Exp) == i) |>
        dplyr::select(Rep, Trait) |>
        filter(!is.na(Trait) & !duplicated(get(Rep))) |>
        dplyr::select(Trait) |>
        unlist() %>%
        sum(.) %>%
        is.na(.)
      ) == F & unique(data1[[Rep]]) > 1) {
        Model2 <-
          as.formula(TraitT ~ as.factor(get(Trat)) + as.factor(get(Rep)))
      } else {
        Model2 <- as.formula(TraitT ~ as.factor(get(Trat)))
      }
      DISCREPANT <- list()
      
      MOD2 <- list()
      for (i in unique(unlist(data2[[Exp]]))) {
        # mod1 = Com NA nos resíduos
        mod1 <-
          lm(Model2,
             data = data2 %>%
               filter(get(Exp) == i),
             na.action = na.omit)
        
        # mod2 = sem NA nos resíduos
        mod2 <-
          lm(Model2,
             data = data2 %>%
               filter(get(Exp) == i),
             na.action = na.exclude)
        MOD2[[i]] <- list(
          mod1 = mod1,
          mod2 = mod2,
          DATA = data2 %>%
            filter(get(Exp) == i)
        )
      }
      
      for (i in unique(unlist(data2[[Exp]]))) {
        # Extraindo o resíduo estudentizado
        MOD2[[i]]$DATA[, paste0("residuals.", Trait)] <-
          rstudent(MOD2[[i]]$mod2, type = "pearson")
        # Extraindo os valores preditos
        MOD2[[i]]$DATA[, paste0("predict.", Trait)] <-
          predict(MOD2[[i]]$mod2)
        # Obtendo os resíduos estudentizados
        # Sem NA's
        MOD2[[i]]$rs <- rstudent(MOD2[[i]]$mod1, type = "pearson")
        
        # Com NA's
        MOD2[[i]]$rs.n.na <-
          rstudent(MOD2[[i]]$mod2, type = "pearson")
      }
      for (i in unique(unlist(data2[[Exp]]))) {
        discrepant <-
          MOD2[[i]]$DATA[which(MOD2[[i]]$DATA[, paste0("residuals.", Trait)] < -3 |
                                 MOD2[[i]]$DATA[, paste0("residuals.", Trait)] > 3), colnames(MOD2[[i]]$DATA)]
        
        DISCREPANT[[i]] <- discrepant
        MOD2[[i]]$DATA[rownames(DISCREPANT[[i]]), "TraitT"] <- NA
        # MOD[[i]]$DATA[rownames(DISCREPANT[[i]]), Trait] <- NA
      }
      
      DATA2 <- list()
      for (i in unique(unlist(data2[[Exp]]))) {
        DATA2[[i]] <- MOD2[[i]]$DATA
      }
      for (i in unique(unlist(data2[[Exp]]))) {
        DATA2[[i]][[Trait]] <- ifelse(is.na(DATA2[[i]]$TraitT), NA,
                                      DATA2[[i]][[Trait]])
        
      }
      data2 <- bind_rows(DATA2)
      
      disc <- disc + 1
      
      if (disc > nDiag) {
        break
      }
    }
    
    
    if (plot_diag2 == TRUE) {
      YP2 <- list()
      LMAX_RS2 <- list()
      LMIN_RS2 <- list()
      LMAX_YP2 <- list()
      LMIN_YP2 <- list()
      
      if (max(MOD[[i]]$rs, na.rm = T) > 3 |
          min(MOD[[i]]$rs, na.rm = T) < -3) {
        DI2 <- "\nAfter removal of discrepant data"
      } else {
        DI2 <- ""
        
      }
      
      for (i in unique(unlist(data2[[Exp]]))) {
        # Histograma
        hist(
          MOD2[[i]]$rs.n.na,
          main = "",
          xlab = "Studentized Residuals",
          ylab = paste("Probability density"),
          freq = F
        )
        title(
          main = paste("Histogram for", Trait, " trait:", i, DI2),
          cex.main = 1,
          sub = SUBt[[i]]
        )
        xm <-
          seq(min(MOD[[i]]$rs, na.rm = T), max(MOD[[i]]$rs, na.rm = T), length = 40)
        ym <- dnorm(xm)
        lines(xm, ym)
        
        # qqNorm
        ## Plotando os pontos
        qqnorm(MOD2[[i]]$rs.n.na,
               main = "")
        # Plotando a reta
        qqline(MOD2[[i]]$rs.n.na,
               main = "",
               xlab = "",
               ylab = "")
        title(
          main = paste("Normality Graph - qqNorm:", Trait, "-", i, DI2),
          cex.main = 1,
          sub = SUBt[[i]]
        )
        # qqPlot
        qqPlot(MOD[[i]]$rs.n.na,
               main = "",
               ylab = "Studentized Residuals")
        title(
          main = paste("Normality Graph - qqPlot:", Trait, "-", i, DI2),
          cex.main = 1,
          sub = SUBt[[i]]
        )
        
        # Checando presença de valores discrepantes
        # Armazenando yp = Valores preditos
        YP2[[i]] <- predict(MOD2[[i]]$mod1, na.rm = T)
        
        # Marcando limites para valores normais
        LMAX_RS2[[i]] <- max(MOD2[[i]]$rs, 3.0, na.rm = T) + 0.1
        LMIN_RS2[[i]] <- min(MOD2[[i]]$rs, -3.0, na.rm = T) - 0.1
        
        # Armazenando os valores preditos
        LMAX_YP2[[i]] <- max(YP2[[i]], na.rm = T) + 0.1
        LMIN_YP2[[i]] <- min(YP2[[i]], na.rm = T) - 0.1
        
        # Plotando o gráfico
        ## Configurando os eixos e títulos
        plot(
          c(LMIN_YP2[[i]], LMAX_YP2[[i]]),
          c(LMIN_RS2[[i]], LMAX_RS2[[i]]),
          "n",
          main = "",
          xlab = expression(paste(hat(y)[predito])),
          ylab = "Studentized Residuals"
        )
        title(
          main = paste("Residual Analysis:", Trait, "-", i, DI2),
          cex.main = 1,
          sub = SUBt[[i]]
        )
        
        ## Plotando as linhas
        abline(h = 0, col = "black")
        abline(h = 3.0, col = "red")
        abline(h = -3.0, col = "red")
        
        ## Plotando os pontos
        points(YP2[[i]], MOD2[[i]]$rs)
      }
    } else {
      DI2 <- ""
      message("plot_diag2 = FALSE: omitting graphs after discrepant analysis")
    }
    
    NORMAL2 <- list()
    ASSIMETRIA2 <- list()
    CURTOSE2 <- list()
    HOMOCED2 <- list()
    # Normalidade e homocedasticidade
    for (i in unique(unlist(data2[[Exp]]))) {
      # Teste de normalidade
      Normal <-
        list(ad.test(MOD2[[i]]$rs), lillie.test(MOD2[[i]]$rs))
      NORMAL2[[i]] <- Normal
      # Assimetria
      Assimetria <- skewness(MOD2[[i]]$rs.n.na, na.rm = T)
      ASSIMETRIA2[[i]] <- Assimetria
      # Curtose
      Curtose <- kurtosis(MOD2[[i]]$rs.n.na, na.rm = T)
      CURTOSE2[[i]] <- Curtose
      
      # Homocedasticidade
      homoced <- leveneTest(TraitT ~ as.factor(get(Trat)),
                            data = data2)
      HOMOCED2[[i]] <- homoced
      
      if (T_BoxCox == TRUE) {
        label(NORMAL2[[i]]) <-
          paste("\nNormality test for ",
                i,
                " after BoxCox transformation",
                DI2)
        label(ASSIMETRIA2[[i]]) <-
          paste("\nSkewness in the data of ",
                i,
                " after BoxCox transformation",
                DI2)
        label(CURTOSE2[[i]]) <-
          paste("\nKurtosis in the data of ",
                i,
                " após transformação BoxCox",
                DI2)
        label(HOMOCED2[[i]]) <-
          paste("\nHomoscedasticity Test for",
                i,
                " after BoxCox transformation",
                DI2)
      } else {
        label(NORMAL2[[i]]) <-
          paste0("\nNormality test for ", i, DI2)
        label(ASSIMETRIA2[[i]]) <-
          paste0("\nSkewness in the data of ", i, DI2)
        label(CURTOSE2[[i]]) <-
          paste0("\nKurtosis in the data of ", i, DI2)
        label(HOMOCED2[[i]]) <-
          paste0("\nHomocedasticity Test for ", i, DI2)
      }
    }
    
    Discrepantes <- suppressMessages(anti_join(data1, data2))
    
    
    if (verbose == TRUE) {
      print(
        list(
          Outliers = Discrepantes,
          `Normality test before discrepant analylis` = NORMAL1,
          `Normality test after discrepant analysis` = NORMAL2,
          `Skewness before discrepant analysis` = ASSIMETRIA1,
          `Skewnwss after discrepant analysis` = ASSIMETRIA2,
          `Kurtosis before discrepant analysis` = CURTOSE1,
          `Kurtosis after discrepant analysis` = CURTOSE2,
          `Homocedasticity test before discrepant analysis` = HOMOCED1,
          `Homocedasticity test after discrepant analysis` = HOMOCED2
        )
      )
    }
    
    Output <-
      list(
        DataOriginals = data1,
        DataAfterDiag = data2[, c(ColumnNames_To_Return)],
        Discrepants = Discrepantes,
        `Normality and Homocedasticity` = list(
          `Normality before diagnostic analysis` = NORMAL1,
          `Normality after diagnostic analysis` = NORMAL2,
          `Skewness before diagnostic analysis` = ASSIMETRIA1,
          `Skewness after diagnostic analysis` = ASSIMETRIA2,
          `Kurtosis before diagnostic analysis` = CURTOSE1,
          `Kurtosis after diagnostic analysis` = CURTOSE2,
          `Homocedasticity before diagnostic analysis` = HOMOCED1,
          `Homocedasticity after diagnostic analysis` = HOMOCED2
        ),
        Lambda = LAMBDA
      )
    return(Output)
  }



#---------------------------- Function for model comparison --------------------------------------#
## 5.1. AIC and -2LogLik values
COMP_MODEL <- function(Traits = c("DAP", "ALT", "VOL"),
                       #Site_Column = "Local",
                       #Age_Column = "Idade",
                       Age = NULL,
                       # Models
                       standard_model = NULL,
                       spatial_model = NULL,
                       competition_model = NULL,
                       competition_spatial_model = NULL) {
  if (is.null(Age)) {
    standard_model2 <- list()
    spatial_model2 <- list()
    competition_model2 <- list()
    competition_spatial_model2 <- list()
    for (i in names(standard_model)) {
      age = "1"
      # Creating lists
      
      # Filling lists
      standard_model2[[i]][[age]] <- standard_model[[i]]
      spatial_model2[[i]][[age]] <- spatial_model[[i]]
      competition_model2[[i]][[age]] <- competition_model[[i]]
      competition_spatial_model2[[i]][[age]] <-
        competition_spatial_model[[i]]
    }
    # Replace lists
    standard_model <- standard_model2
    spatial_model <- spatial_model2
    competition_model <- competition_model2
    competition_spatial_model <- competition_spatial_model2
    # Remove lists
    rm(
      standard_model2,
      spatial_model2,
      competition_model2,
      competition_spatial_model2
    )
    Age = age
    
  }
  # Names models
  models <-
    c(
      "standard_model",
      "spatial_model",
      "competition_model",
      "competition_spatial_model"
    )
  # Create list with models
  MODELS <- as.list(mget(models))
  
  # if (class(Model)[1] != "breedR") {
  #       stop("This is not a breedR adjusted model")
  # }
  
  # Remove non-existent models
  TrueModel <- lapply(MODELS, function(x)
    ! is.null(unlist(x)))
  TrueModel <- names(TrueModel[which(TrueModel == T)])
  
  # Extract
  AIC_LIST <- list()
  for (i in TrueModel) {
    for (j in Traits) {
      for (k in Age) {
        AIC_LIST[[i]][[j]][[k]] <- list(
          AIC = MODELS[[i]][[j]][[k]][["fit"]][["AIC"]],
          Loglik = MODELS[[i]][[j]][[k]][["fit"]][["-2logL"]],
          Df = summary(MODELS[[i]][[j]][[k]])$var |>
            unlist() |>
            unique() |>
            length()
        )
      }
    }
  }
  
  DIF_LogLik <- list()
  Df_MODELS <- list()
  for (i in Traits) {
    for (j in Age) {
      # Extract Loglik
      DIF_LogLik[[i]][[j]] <-
        lapply(AIC_LIST, function(x)
          x[[i]][[j]][["Loglik"]]) |>
        unlist(use.names = T) %>%
        outer(., ., FUN = "-") |>
        abs()
      # Extract Degree of freedom
      Df_MODELS <-
        lapply(AIC_LIST, function(x)
          x[[i]][[j]][["Df"]]) |>
        unlist()
    }
  }
  
  for (i in TrueModel) {
    for (j in Traits) {
      for (k in Age) {
        # LRT value
        AIC_LIST[[i]][[j]][[k]][["LRT_value"]] <-
          DIF_LogLik[[j]][[k]]["standard_model", ][i]
        # LRT p.value
        AIC_LIST[[i]][[j]][[k]][["LRT_p.value"]] <-
          pchisq(
            DIF_LogLik[[j]][[k]]["standard_model", ][i],
            df = abs(Df_MODELS["standard_model"] - Df_MODELS[i]),
            lower.tail = F
          )
      }
    }
  }
  
  MODEL_COMPARISON <- list()
  for (i in TrueModel) {
    for (j in Traits) {
      for (k in Age) {
        MODEL_COMPARISON[[i]][[j]][[k]] <-
          data.frame(AIC_LIST[[i]][[j]][[k]]) |>
          rownames_to_column(var = "Model") |>
          mutate(Traits = j,
                 Age = k)
      }
    }
  }
  
  #rm(M_LAPPLY)
  M_LAPPLY <- lapply(MODEL_COMPARISON, function(x)
    lapply(x, function(y)
      bind_rows(y))) |>
    lapply(function(x)
      bind_rows(x)) |>
    bind_rows()
  
  return(M_LAPPLY)
}

#--------------------------- Function for extract heritability and correlations ------------------#

Extract_h2a <- function(Model = NULL,
                        random_effect = NULL,
                        model_type = "std") {
  Var = summary(Model)$var
  # For ar1Comp Model
  if (model_type == "ar1Comp") {
    test_model <- names(Var)
    if (is.null(test_model)) {
      stop("model_type is not 'ar1Comp' model, please, set properly!")
    }
    if (!"spatial" %in% names(Var)) {
      stop("model_type is not 'ar1Comp' model, please, set properly!")
    }
    if (!is.null(random_effect)) {
      Names_Var <- c(names(Var), "competition", "corr_gen_comp")
      Variances_Comp <- data.frame(
        Stats = Names_Var,
        Estimates = c(
          Var[[random_effect]],
          Var$genetic[1, 1],
          Var$pec,
          Var$spatial,
          Var$Residual,
          Var$genetic[2, 2],
          Var$genetic[1, 2]
        ),
        S.E = as.numeric(NA)
      ) |>
        `rownames<-`(Names_Var)
      # Correlation between direct and indirect genetic effect
      Corr_gen_comp <-
        Variances_Comp["corr_gen_comp", "Estimates"] / (sqrt(Variances_Comp["genetic", "Estimates"]) * sqrt(Variances_Comp["competition", "Estimates"]))
      Variances_Comp["corr_gen_comp", "Estimates"] <- Corr_gen_comp
      # spatial variance
      Spatial <- Variances_Comp["spatial", "Estimates"]
      
      # c2_random_effect
      C2 <- list()
      for (i in random_effect) {
        C2[[i]] <- data.frame(
          Stats = paste0("c2_", i),
          Estimates = Variances_Comp[random_effect, "Estimates"] / (sum(abs(
            Variances_Comp$Estimates
          )) - abs(Corr_gen_comp) - Spatial),
          S.E = NA
        )
        
      }
      c2_random_effect <- bind_rows(C2)
      
      # Variance components
      # h2a
      h2a <-
        c("h2a", Variances_Comp["genetic", "Estimates"] / (sum(abs(
          Variances_Comp$Estimates
        )) - abs(Corr_gen_comp) - Spatial), NA)
      Data2 <- rbind(Variances_Comp, c2_random_effect, h2a) |>
        `rownames<-`(NULL)
      
    } else {
      Names_Var <- c(names(Var), "competition", "corr_gen_comp")
      Variances_Comp <- data.frame(
        Stats = Names_Var,
        Estimates = c(
          Var$genetic[1, 1],
          Var$pec,
          Var$spatial,
          Var$Residual,
          Var$genetic[2, 2],
          Var$genetic[1, 2]
        ),
        S.E = as.numeric(NA)
      ) |>
        `rownames<-`(Names_Var)
      # Correlation between direct and indirect genetic effect
      Corr_gen_comp <-
        Variances_Comp["corr_gen_comp", "Estimates"] / (sqrt(Variances_Comp["genetic", "Estimates"]) * sqrt(Variances_Comp["competition", "Estimates"]))
      Variances_Comp["corr_gen_comp", "Estimates"] <- Corr_gen_comp
      # spatial variance
      Spatial <- Variances_Comp["spatial", "Estimates"]
      
    }
    # Variance components
    # h2a
    h2a <-
      c("h2a", Variances_Comp["genetic", "Estimates"] / (sum(abs(
        Variances_Comp$Estimates
      )) - abs(Corr_gen_comp) - Spatial), NA)
    
    Data2 <- rbind(Variances_Comp, h2a) |>
      `rownames<-`(NULL)
  } else {
    if (model_type == "Comp") {
      if (is.null(names(Var))) {
        stop("model_type is not 'Comp' model, please, set properly!")
      }
      if ("spatial" %in% names(Var)) {
        stop("model_type is not 'Comp' model, please, set properly!")
      }
      if (!is.null(random_effect)) {
        Names_Var <- c(names(Var), "competition", "corr_gen_comp")
        Variances_Comp <- data.frame(
          Stats = Names_Var,
          Estimates = c(
            Var[[random_effect]],
            Var$genetic[1, 1],
            Var$pec,
            Var$Residual,
            Var$genetic[2, 2],
            Var$genetic[1, 2]
          ),
          S.E = as.numeric(NA)
        ) |>
          `rownames<-`(Names_Var)
        # Correlation between direct and indirect genetic effect
        Corr_gen_comp <-
          Variances_Comp["corr_gen_comp", "Estimates"] / (sqrt(Variances_Comp["genetic", "Estimates"]) * sqrt(Variances_Comp["competition", "Estimates"]))
        Variances_Comp["corr_gen_comp", "Estimates"] <-
          Corr_gen_comp
        # c2_random_effect
        C2 <- list()
        for (i in random_effect) {
          C2[[i]] <- data.frame(
            Stats = paste0("c2_", i),
            Estimates = Variances_Comp[random_effect, "Estimates"] / (sum(Variances_Comp$Estimates) - Corr_gen_comp),
            S.E = NA
          )
          
        }
        c2_random_effect <- bind_rows(C2)
        # h2a
        h2a <-
          c("h2a", Variances_Comp["genetic", "Estimates"] / (sum(Variances_Comp$Estimates) - Corr_gen_comp), NA)
        
        Data2 <- rbind(Variances_Comp, c2_random_effect, h2a) |>
          `rownames<-`(NULL)
        
      } else{
        Names_Var <- c(names(Var), "competition", "corr_gen_comp")
        Variances_Comp <- data.frame(
          Stats = Names_Var,
          Estimates = c(
            Var$genetic[1, 1],
            Var$pec,
            Var$Residual,
            Var$genetic[2, 2],
            Var$genetic[1, 2]
          ),
          S.E = as.numeric(NA)
        ) |>
          `rownames<-`(Names_Var)
        # Correlation between direct and indirect genetic effect
        Corr_gen_comp <-
          Variances_Comp["corr_gen_comp", "Estimates"] / (sqrt(Variances_Comp["genetic", "Estimates"]) * sqrt(Variances_Comp["competition", "Estimates"]))
        Variances_Comp["corr_gen_comp", "Estimates"] <-
          Corr_gen_comp
        # c2_random_effect
        C2 <- list()
        for (i in random_effect) {
          C2[[i]] <- data.frame(
            Stats = paste0("c2_", i),
            Estimates = Variances_Comp[i, "Estimates"] / (sum(Variances_Comp$Estimates) - Corr_gen_comp),
            S.E = NA
          )
          
        }
        c2_random_effect <- bind_rows(C2)
        
        # h2a
        h2a <-
          c("h2a", Variances_Comp["genetic", "Estimates"] / (sum(Variances_Comp$Estimates) - Corr_gen_comp), NA)
        
        Data2 <- rbind(Variances_Comp, h2a) |>
          `rownames<-`(NULL)
      }
    } else {
      # AR1 Model
      if (model_type == "ar1") {
        if (is.null(dimnames(Var))) {
          stop("model_type is not 'ar1' model, please, set properly!")
        }
        if (!"spatial" %in% dimnames(Var)[[1]]) {
          stop("model_type is not 'ar1' model, please, set properly!")
        }
        if (!is.null(random_effect)) {
          Names_Var <- c(dimnames(Var)[[1]], "rho_r", "rho_c")
          RHO <-
            summary(Model)$rho[which.max(summary(Model)$rho$loglik), ]
          Variances_Comp <- data.frame(
            Stats = Names_Var,
            Estimates = c(Var[random_effect, "Estimated variances"],
                          Var["genetic", "Estimated variances"],
                          Var["spatial", "Estimated variances"],
                          Var["Residual", "Estimated variances"],
                          RHO[["rho_r"]],
                          RHO[["rho_c"]]),
            S.E = as.numeric(NA)
          ) |>
            `rownames<-`(Names_Var)
          
          # spatial variance
          Spatial <- Variances_Comp["spatial", "Estimates"]
          
          # c2_random_effect
          C2 <- list()
          for (i in random_effect) {
            C2[[i]] <- data.frame(
              Stats = paste0("c2_", i),
              Estimates = Variances_Comp[i, "Estimates"] / (sum(
                abs(Variances_Comp$Estimates)
              ) - Spatial),
              S.E = NA
            )
            
          }
          c2_random_effect <- bind_rows(C2)
          
          # Variance components
          # h2a
          h2a <-
            c("h2a", Variances_Comp["genetic", "Estimates"] / (sum(
              abs(Variances_Comp$Estimates)
            ) - Spatial), NA)
          Data2 <- rbind(Variances_Comp, c2_random_effect, h2a) |>
            `rownames<-`(NULL)
          
        } else {
          Names_Var <- c(dimnames(Var)[[1]], "rho_r", "rho_c")
          RHO <-
            summary(Model)$rho[which.max(summary(Model)$rho$loglik), ]
          Variances_Comp <- data.frame(
            Stats = Names_Var,
            Estimates = c(Var["genetic", "Estimated variances"],
                          Var["spatial", "Estimated variances"],
                          Var["Residual", "Estimated variances"],
                          RHO[["rho_r"]],
                          RHO[["rho_c"]]),
            S.E = as.numeric(NA)
          ) |>
            `rownames<-`(Names_Var)
          
          # spatial variance
          Spatial <- Variances_Comp["spatial", "Estimates"]
          
          # Variance components
          # h2a
          h2a <-
            c("h2a", Variances_Comp["genetic", "Estimates"] / (sum(
              abs(Variances_Comp$Estimates)
            ) - Spatial), NA)
          Data2 <- rbind(Variances_Comp, h2a) |>
            `rownames<-`(NULL)
        }
      } else{
        if (is.null(dimnames(Var))) {
          stop("model_type is not 'std' model, please, set properly!")
        }
        if ("spatial" %in% dimnames(Var)[[1]]) {
          stop("model_type is not 'std' model, please, set properly!")
        }
        # Standard model
        if (!is.null(random_effect)) {
          Names_Var <- dimnames(Var)[[1]]
          Variances_Comp <- data.frame(
            Stats = Names_Var,
            Estimates = c(Var[random_effect, "Estimated variances"],
                          Var["genetic", "Estimated variances"],
                          Var["Residual", "Estimated variances"]),
            S.E = as.numeric(NA)
          ) |>
            `rownames<-`(Names_Var)
          
          # c2_random_effect
          C2 <- list()
          for (i in random_effect) {
            C2[[i]] <- data.frame(
              Stats = paste0("c2_", i),
              Estimates = Variances_Comp[i, "Estimates"] / (sum(
                abs(Variances_Comp$Estimates)
              )),
              S.E = NA
            )
            
          }
          c2_random_effect <- bind_rows(C2)
          
          
          # Variance components
          # h2a
          h2a <-
            c("h2a", Variances_Comp["genetic", "Estimates"] / (sum(
              abs(Variances_Comp$Estimates)
            )), NA)
          Data2 <- rbind(Variances_Comp, c2_random_effect, h2a) |>
            `rownames<-`(NULL)
          
        } else {
          Names_Var <- dimnames(Var)[[1]]
          Variances_Comp <- data.frame(
            Stats = Names_Var,
            Estimates = c(Var["genetic", "Estimated variances"],
                          Var["Residual", "Estimated variances"]),
            S.E = as.numeric(NA)
          ) |>
            `rownames<-`(Names_Var)
          
          # Variance components
          # h2a
          h2a <-
            c("h2a", Variances_Comp["genetic", "Estimates"] / (sum(
              abs(Variances_Comp$Estimates)
            )), NA)
          Data2 <- rbind(Variances_Comp, h2a) |>
            `rownames<-`(NULL)
        }
      }
    }
  }
  return(Data2)
}

#------------ Function for perform deviance analysis from remlf90 function -----#

Deviance_BreedR <-
  function(Trait = Trait,
           Model = breedR_model,
           Data = Progenies,
           Pedigree = NULL,
           Method = "EM",
           model_type = "std_animal") {
    if (sum(stringr::str_count(names(Model[["mf"]]), Trait)) < 1 &
        sum(stringr::str_count(names(Model[["mf"]]), "get")) < 1) {
      stop(paste(
        Trait,
        "trait not found. Please check the evaluated trait in the model"
      ))
    }
    if (!model_type %in% c("random", "std_animal")) {
      stop("model_type is not appropriate. Please, choose between random or std_animal model")
    }
    
    # Detect fixed effects
    fixed_effect <- names(Model[["fixed"]])
    # Detect random effects
    random_effect <- names(Model[["ranef"]])
    # random effect without genetic effect
    random_effect2 <-
      random_effect[str_detect(random_effect, "genetic", negate = T)]
    
    # Function to transform variables in factors
    to.factors <- function(df, variables) {
      for (variable in variables) {
        df[[variable]] <- as.factor(df[[variable]])
      }
      return(df)
    }
    
    # Convert effects in the model in factors
    Data <- to.factors(Data, c(fixed_effect, random_effect2))
    
    if (is.null(Pedigree)) {
      ped <<- get(Model[["call"]][["genetic"]][["pedigree"]])
    } else {
      ped <<- Pedigree
    }
    
    
    DEV <- list()
    MOD <- list()
    #FIT.2logL <- list()
    
    for (i in random_effect2) {
      # loop for model_type
      # if model_type == 'random':
      if (model_type == "random") {
        # Update models
        Model2 <- remlf90(
          as.formula(Model[["call"]][["fixed"]]),
          random = update.formula(Model[["call"]][["random"]],
                                  formula(paste("~ . ", "-", i))),
          dat = Data,
          method = Method
        )
        
      } else {
        # if model_type == 'animal':
        if (sum(stringr::str_count(names(Model[["mf"]]), Trait)) < 1 &
            sum(stringr::str_count(names(Model[["mf"]]), "get")) > 0) {
          model_fix <-  paste(Trait, paste("~", sub("get(+$)", "", Model[["call"]][["fixed"]])[3]), collapse = "")
        } else {
          model_fix <- Model[["call"]][["fixed"]]
        }
        
        
        Model2 <- remlf90(
          as.formula(model_fix),
          random = update.formula(Model[["call"]][["random"]],
                                  formula(paste("~ . ", "-", i))),
          genetic = list(
            model = "add_animal",
            pedigree = ped,
            id = Model[["call"]][["genetic"]][["id"]]
          ),
          dat = Data,
          method = Method
        )
        # Model without genetic effect
        
        Model3 <- remlf90(
          as.formula(model_fix),
          random = as.formula(Model[["call"]][["random"]]),
          dat = Data,
          method = Method
        )
      }
      # Extract -2logL
      TestStat <-
        Model2[["fit"]][["-2logL"]] - Model[["fit"]][["-2logL"]]
      # Extract parameters differences between models
      DifGL <-
        length(names(Model[["ranef"]])) - length(names(Model2[["ranef"]]))
      # Obtain the p-value
      p.value <- pchisq(TestStat, df = DifGL, lower.tail = FALSE)
      # Building the deviance table
      DataDev <-
        data.frame(Efeito = i,
                   LRT = TestStat,
                   p.valor = p.value)
      # FIT.2logL[[i]] <- Model2[["fit"]][["-2logL"]]
      DEV[[i]] <- DataDev
      
    }
    Deviance <- bind_rows(DEV)
    
    DevGenetic <-
      c(
        "genetic",
        Model3[["fit"]][["-2logL"]] - Model[["fit"]][["-2logL"]],
        pchisq(
          Model3[["fit"]][["-2logL"]] - Model[["fit"]][["-2logL"]],
          df = length(names(Model[["ranef"]])) - length(names(Model3[["ranef"]])),
          lower.tail = FALSE
        )
      )
    #DF <- length(names(Model[["ranef"]])) - length(names(Model3[["ranef"]]))
    Deviance[3, ] <- DevGenetic
    return(Deviance)
  }



#--------------------------- Function for get breeding values ------------------#

BV_BreedR <- function(data = data,
                      Model = NULL,
                      NamesID = "ID2",
                      Pedigree = Pedigree,
                      dominanceName = NULL,
                      ColumnFamilies = "Prog",
                      Rank = T,
                      NamesReturnDataset = c("ID2", "Mae", "Bloco", "Arvore", "DAP", "ALT"),
                      model_type = "std",
                      random_effect = NULL,
                      fixed_effect_av = "Intercept") {
  if (class(Model)[1] != "breedR") {
    stop("This is not a breedR adjusted model")
  }
  
  # Map original individuals' ID
  MapPed <- attr(Pedigree, "map")
  
  # Get the original individuals' ID
  labelPed <- match(get_pedigree(Model)@label, MapPed)
  nParents <-
    as.matrix(summary(as.data.frame(get_pedigree(Model))$dam))["Max.", 1]
  nId <-
    as.matrix(summary(as.data.frame(get_pedigree(Model))$self))["Max.", 1] - nParents
  # Average
  mu <- mean(fixef(Model)[[fixed_effect_av]])
  
  # Standard and AR1 model
  if (model_type == "std" | model_type == "ar1") {
    if (model_type == "std" & !"genetic" %in% names(ranef(Model))) {
      stop("model_type is not 'std' model, please, set properly!")
    }
    if (model_type == "ar1" & !"genetic" %in% names(ranef(Model))) {
      stop("model_type is not 'ar1' model, please, set properly!")
    }
    # BV Parents
    BV_Parents <-
      Model[["ranef"]][["genetic"]][[1]][1:nParents, ] |>
      `rownames<-`(labelPed[1:nParents]) |>
      mutate(`u+a` = mu + value) |>
      rownames_to_column(var = "Family") |>
      `colnames<-`(c("Family", "a", "s.e", "u+a"))
    # BV Individuals
    BV_Ind <-
      Model[["ranef"]][["genetic"]][[1]][nParents + 1:nId, ] |>
      `rownames<-`(labelPed[nParents + 1:nId]) |>
      mutate(`u+a` = mu + value) |>
      rownames_to_column(var = NamesID) |>
      `colnames<-`(c(NamesID, "a", "s.e", "u+a"))
    # Genetic variance
    Va <- summary(Model)$var["genetic", "Estimated variances"]
    
    if (is.null(random_effect)) {
      BV_random_effect = NULL
    } else {
      # BV random effect
      BV_random_effect <- Model[["ranef"]][[random_effect]][[1]]
    }
  } else {
    # Comp and AR1-Comp model
    if (model_type == "Comp" | model_type == "ar1Comp") {
      if (model_type == "Comp" &
          !"genetic_direct" %in% names(ranef(Model))) {
        stop("model_type is not 'Comp' model, please, set properly!")
      }
      if (model_type == "ar1" &
          !"genetic_direct" %in% names(ranef(Model))) {
        stop("model_type is not 'ar1Comp' model, please, set properly!")
      }
      # Direct effect
      # BV Parents
      BV_Parents <-
        Model[["ranef"]][["genetic_direct"]][[1]][1:nParents, ] |>
        `rownames<-`(labelPed[1:nParents]) |>
        mutate(`u+a` = mu + value) |>
        rownames_to_column(var = "Family") |>
        `colnames<-`(c("Family", "a", "s.e", "u+a"))
      # BV Individuals
      BV_Ind <-
        Model[["ranef"]][["genetic_direct"]][[1]][nParents + 1:nId, ] |>
        `rownames<-`(labelPed[nParents + 1:nId]) |>
        mutate(`u+a` = mu + value) |>
        rownames_to_column(var = NamesID) |>
        `colnames<-`(c(NamesID, "a", "s.e", "u+a"))
      
      # Competition effect
      BV_Competition_Parents <-
        Model[["ranef"]][["genetic_competition"]][[1]][1:nParents, ] |>
        `rownames<-`(labelPed[1:nParents]) |>
        mutate(`u+a` = mu + value) |>
        rownames_to_column(var = "Family") |>
        `colnames<-`(c(
          "Family",
          "competition",
          "s.e_competition",
          "u+competition"
        ))
      # BV Individuals
      BV_Competition_Ind <-
        Model[["ranef"]][["genetic_competition"]][[1]][nParents + 1:nId, ] |>
        `rownames<-`(labelPed[nParents + 1:nId]) |>
        mutate(`u+a` = mu + value) |>
        rownames_to_column(var = NamesID) |>
        `colnames<-`(c(
          NamesID,
          "competition",
          "s.e_competition",
          "u+competition"
        ))
      
      # Genetic variance
      Va <-
        summary(Model)$var$genetic["genetic_direct", "genetic_direct"]
      if (is.null(random_effect)) {
        BV_random_effect = NULL
      } else {
        # BV random effect
        BV_random_effect <- Model[["ranef"]][[random_effect]][[1]]
      }
    }
  }
  # Accuracy from family model
  Accuracy.Family = mean(sqrt(1 - BV_Parents$s.e ^ 2 / Va))
  # Accuracy from ID model
  Accuracy.Ind <- mean(sqrt(1 - BV_Ind$s.e ^ 2 / Va), na.rm = T)
  # Change BV_Ind class
  class(BV_Ind[[NamesID]]) <- class(data[, NamesID])
  Data_Total <- left_join(data[, NamesReturnDataset],
                          BV_Ind, by = NamesID)
  
  # attr(BV_random_effect, "class") <- c("data.frame", "BV_BreedR")
  # attr(BV_Parents, "class") <- c("data.frame", "BV_BreedR")
  # attr(BV_Ind, "class") <- c("data.frame", "BV_BreedR")
  # attr(Data_Total, "class") <- c("data.frame", "BV_BreedR")
  #
  Results <- list(
    BV_random_effect = BV_random_effect,
    BV_Parents = BV_Parents,
    BV_Ind = BV_Ind,
    Data_Total = Data_Total,
    Accuracy.Family = Accuracy.Family,
    Accuracy.Ind = Accuracy.Ind
  )
  return(Results)
}

#--------------------------- Function for extract estimates from lmer model ---------------#
ExtractLmer <- function(Model = Model,
                        Trat = "tratamento",
                        Rep = "bloco",
                        nRep = 6,
                        nPlant = 6,
                        Plot = TRUE,
                        random = NULL) {
  if (class(Model)[1] != "lmerModLmerTest") {
    stop("This is not a lmerModLmerTest adjusted model")
  }
  
  # Extract variance components
  VarComp <- as.data.frame(VarCorr(Model)) %>%
    dplyr::select(!c("var1", "var2")) %>%
    `colnames<-`(c("FV", "Variance", "Standard Deviation")) %>%
    `rownames<-`(., .$FV)
  
  Plot <- ifelse(Plot == FALSE, NA, Plot)
  
  EST <- tibble(
    Vg = VarComp[Trat, "Variance"],
    Vplot = ifelse(Plot == TRUE, VarComp[paste0(Trat, sep = ":", Rep), "Variance"], Plot),
    Vrandom = ifelse(is.null(random), NA, VarComp[random, "Variance"]),
    Vres = VarComp["Residual", "Variance"],
    Vpheno = sum(Vg, Vplot, Vrandom, Vres, na.rm = T),
    `h2g` = Vg / Vpheno,
    `h2m` = Vg / sum(Vg, (Vplot / nRep), (Vres / (nRep * nPlant)), na.rm = T),
    `c2p` = ifelse(is.na(Vplot), NA, Vplot / Vpheno)
  ) %>%
    t() %>%
    as.data.frame() %>%
    `colnames<-`("Estimates") %>%
    rownames_to_column(var = "Components")
  return(EST)
}

#----------------- Function for extract BLUPs from lmer ------------#
BV_lmer <- function(Model = Model,
                    Trat = "tratamento",
                    FixEffect = "bloco",
                    Rank = TRUE) {
  if (class(Model)[1] != "lmerModLmerTest") {
    stop("This is not a lmerModLmerTest adjusted model")
  }
  
  # Estimated Marginal Mean
  mu <- emmeans::emmeans(Model, specs = FixEffect) |>
    data.frame() |>
    summarise(mu = mean(emmean)) |>
    unlist()
  
  # BLUPs
  blups <- ranef(Model)[[Trat]] |>
    `colnames<-`("a") |>
    rownames_to_column(Trat) |>
    mutate(`u+a` = mu + a)
  
  # Arrange?
  if (Rank == TRUE) {
    blups <- blups |>
      arrange(desc(a))
  } else {
    blups <- blups
  }
  return(blups)
}

#---------------------- Function for thinning strategies --------------------------
Thinning_BreedR <- function(BV_Column = "a_total",
                            Trait = Trait,
                            BV_fam = BV_Family,
                            Data_Total = BV$Data_Total,
                            Family_Data_Total = Family,
                            Bloc_Column = "Block",
                            nGroups = 4,
                            label.group.y = c(1, 1, 1, 1),
                            Plot.Rank = TRUE,
                            save_plot_rank = TRUE,
                            IS = NULL,
                            id = "ID",
                            save_table_xlsx = TRUE) {
  # if (exists("nGroups", mode = "any")) {
  #   nGroups = readline(prompt = "Enter with the number of groups for thinning strategies:")
  # }
  
  if (!nGroups %in% 2:4) {
    stop("nGroups should be: 2, 3 or 4")
  }
  
  BV_fam <- BV_fam |>
    arrange(desc(get(BV_Column))) |>
    mutate(Family = factor(Family, levels = Family))
  
  
  
  # Families with positive BV
  MajorBV_fam <- BV_fam |>
    filter(get(BV_Column) >= 0)
  
  # Familie with BV near of zero
  Fam_Zero <- which.min(MajorBV_fam[[BV_Column]])
  # Families with negative BV
  MinorBV_fam <- BV_fam |>
    filter(get(BV_Column) < 0)
  
  # Found inflexion point in the best families
  Inflex1 <- RootsExtremaInflections::inflexi(
    x = 1:nrow(MajorBV_fam),
    y = MajorBV_fam[[BV_Column]],
    nt = 2,
    i1 = 1,
    i2 = nrow(MajorBV_fam),
    plots = F,
  )
  
  Inflexi_major <- Inflex1$finfl[1]
  
  # Found inflexion point in the worst families
  Inflex2 <- RootsExtremaInflections::inflexi(
    x = 1:nrow(MinorBV_fam),
    y = MinorBV_fam[[BV_Column]],
    nt = 2,
    i1 = 1,
    i2 = nrow(MinorBV_fam),
    plots = F,
  )
  
  Inflexi_minor <- Fam_Zero + Inflex2$finfl[1]
  Col1 = "gray70"
  Col2 = "blue"
  Col3 = "red"
  # line in Fam_Zero
  nGroup2 <-
    geom_vline(
      xintercept = Fam_Zero,
      linetype = 4,
      colour = Col1,
      size = 1
    )
  # nGroup = 3
  nGroup3 <-
    list(
      geom_vline(
        xintercept = Fam_Zero,
        linetype = 4,
        colour = Col1,
        size = 1
      ),
      geom_vline(
        xintercept = Inflexi_major,
        linetype = 4,
        colour = Col2,
        size = 1
      )
    )
  # nGroup = 4
  nGroup4 <-
    list(
      geom_vline(
        xintercept = Fam_Zero,
        linetype = 4,
        colour = Col1,
        size = 1
      ),
      geom_vline(
        xintercept = Inflexi_major,
        linetype = 4,
        colour = Col2,
        size = 1
      ),
      geom_vline(
        xintercept = Inflexi_minor,
        linetype = 4,
        colour = Col3,
        size = 1
      )
    )
  
  # if 2 groups:
  if (nGroups ==  2) {
    nGroupFinal <- nGroup2
    AnnotateGroup <- list(
      annotate(
        "text",
        x = Fam_Zero / 2,
        y = label.group.y[1],
        label = "G1"
      ),
      annotate(
        "text",
        x = (Fam_Zero + (nrow(BV_fam) - Fam_Zero) / 2),
        y = label.group.y[2],
        label = "G2"
      )
    )
  } else {
    # if 3 groups
    if (nGroups == 3) {
      nGroupFinal <- nGroup3
      AnnotateGroup <- list(
        annotate(
          "text",
          x = Inflexi_major / 2,
          y = label.group.y[1],
          label = "G1"
        ),
        annotate(
          "text",
          x = (Fam_Zero - (Fam_Zero - Inflexi_major) / 2),
          y = label.group.y[2],
          label = "G2"
        ),
        annotate(
          "text",
          x = Fam_Zero + Fam_Zero / 2,
          y = label.group.y[3],
          label = "G3"
        )
      )
    } else {
      # if 4 Groups
      nGroupFinal <- nGroup4
      AnnotateGroup <- list(
        annotate(
          "text",
          x = Inflexi_major / 2,
          y = label.group.y[1],
          label = "G1"
        ),
        annotate(
          "text",
          x = (Fam_Zero - (Fam_Zero - Inflexi_major) / 2),
          y = label.group.y[2],
          label = "G2"
        ),
        annotate(
          "text",
          x = Fam_Zero + (Inflexi_minor - Fam_Zero) / 2,
          y = label.group.y[3],
          label = "G3"
        ),
        annotate(
          "text",
          x = Inflexi_minor + (nrow(BV_fam) - Inflexi_minor) / 2,
          y = label.group.y[4],
          label = "G4"
        )
      )
    }
  }
  
  # Plot BV vs Rank from BV
  families_rank <-
    ggplot(BV_fam, aes(x = Family, y = get(BV_Column))) +
    geom_point(size = 2) +
    ggtitle(label = paste("Thinning Strategies based on", Trait, "Trait")) +
    ylab(paste0("Additive genetic effect (", BV_Column, ") from ", Trait, " trait")) + xlab("Families") +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "gray60",
      size = 0.5
    ) +
    nGroupFinal +
    AnnotateGroup
  if (Plot.Rank == TRUE) {
    plot(families_rank)
  }
  
  # Map group for trees in dataset that have BV and experiment information
  Data_Total$Group <- NA
  # G1
  if (nGroups == 2) {
    Data_Total[Data_Total[, Family_Data_Total] %in% as.numeric(as.character(BV_fam[1:Fam_Zero, "Family"])), "Group"] <-
      "G1"
  } else {
    Data_Total[Data_Total[, Family_Data_Total] %in% as.numeric(as.character(BV_fam[1:Inflexi_major, "Family"])), "Group"] <-
      "G1"
  }
  
  # G2
  if (nGroups == 2) {
    Data_Total[Data_Total[, Family_Data_Total] %in% as.numeric(as.character(BV_fam[(Fam_Zero + 1):nrow(BV_fam), "Family"])), "Group"] <-
      "G2"
  } else {
    Data_Total[Data_Total[, Family_Data_Total] %in% as.numeric(as.character(BV_fam[(Inflexi_major + 1):Fam_Zero, "Family"])), "Group"] <-
      "G2"
  }
  
  # G3
  if (nGroups == 3) {
    Data_Total[Data_Total[, Family_Data_Total] %in% as.numeric(as.character(BV_fam[(Fam_Zero + 1):nrow(BV_fam), "Family"])), "Group"] <-
      "G3"
  } else {
    if (nGroups == 4) {
      Data_Total[Data_Total[, Family_Data_Total] %in% as.numeric(as.character(BV_fam[(Fam_Zero + 1):Inflexi_minor, "Family"])), "Group"] <-
        "G3"
    }
    
  }
  
  # G4
  if (nGroups == 4) {
    Data_Total[Data_Total[, Family_Data_Total] %in% as.numeric(as.character(BV_fam[(Inflexi_minor + 1):nrow(BV_fam), "Family"])), "Group"] <-
      "G4"
  }
  
  # Map Plot
  Data_Total <- Data_Total |>
    mutate(Plot = paste0(get(Family_Data_Total), sep = ":", get(Bloc_Column)))
  
  # Combination for thinning strategies
  AllComb <- gtools::combinations(
    n = as.numeric(nGroups) + 1,
    r = as.numeric(nGroups),
    v = (as.numeric(nGroups)):0,
    set = F,
    repeats.allowed = T
  )
  
  AllComb <- AllComb[-nrow(AllComb), ]
  colnames(AllComb) <- paste0("G", 1:nGroups)
  rownames(AllComb) <- seq.int(nrow(AllComb))
  
  
  # Mutate Family_Data_Total
  Data_Total <- Data_Total |>
    mutate(Family = as.character(get(Family_Data_Total)))
  
  
  # BW strategies: Selection between and within families
  gg.BW = list()
  G.BW <- list()
  for (k in sort(unique(Data_Total$Group))) {
    for (i in seq_along(AllComb[, 1])) {
      g.bw <- Data_Total %>%
        filter(Group == k) %>%
        dplyr::group_by(Family, Plot) %>%
        arrange(Family, desc(a))
      if (AllComb[i, which(colnames(AllComb) == k)] == 0) {
        gg.BW[[k]][[i]] <-
          dplyr::slice(g.bw, n = AllComb[i, which(colnames(AllComb) == k)])
      } else {
        gg.BW[[k]][[i]] <-
          dplyr::slice_head(g.bw, n = AllComb[i, which(colnames(AllComb) == k)])
        
      }
      
      if (nGroups == 2) {
        G.BW[[i]] <-
          rbind(gg.BW[["G1"]][[i]], gg.BW[["G2"]][[i]])
      } else {
        if (nGroups == 3) {
          G.BW[[i]] <-
            rbind(gg.BW[["G1"]][[i]], gg.BW[["G2"]][[i]], gg.BW[["G3"]][[i]])
        } else {
          G.BW[[i]] <-
            rbind(gg.BW[["G1"]][[i]], gg.BW[["G2"]][[i]], gg.BW[["G3"]][[i]], gg.BW[["G4"]][[i]])
        }
      }
      
      
      G.BW[[i]][, "BV"] <-
        G.BW[[i]]$a + mean(Data_Total[[Trait]], na.rm = T)
      G.BW[[i]][, "Strategy"] <- "BW"
    }
  }
  # Set names in the strategies
  G.BW <- setNames(G.BW, seq_along(G.BW))
  
  # BI Strategies: Selection of the best individuals regardless of family that he belongs
  ## Second group strategy
  if (is.null(IS)) {
    IC.TOP <- c(seq(0.01, 0.21, 0.02), 0.3, 0.4, 0.5)
    message(
      "The following IS are being used on the BI strategy:\n0.01 0.03 0.05 0.07 0.09 0.11 0.13 0.15 0.17 0.19 0.21 0.30 0.40 0.50"
    )
  } else {
    message(paste("The following IS are being used on the BI strategy:\n", IS))
    IC.TOP <- IS
  }
  
  SeleTop.BI <- list()
  for (i in seq_along(IC.TOP)) {
    seletop <- Data_Total %>%
      group_by(Family, Plot) %>%
      arrange(desc(a)) %>%
      head(n = round(nrow(Data_Total) * IC.TOP[i]))
    seletop$BV <- seletop$a + mean(Data_Total[[Trait]], na.rm = T)
    SeleTop.BI[[i]] <- seletop
    SeleTop.BI[[i]]["Strategy"] <- "BI"
  }
  SeleTop.BI <-
    setNames(SeleTop.BI, seq(1 + length(G.BW), (length(G.BW) + length(SeleTop.BI))))
  
  # Combining BW and BI strategies
  Strategies <- append(G.BW, SeleTop.BI)
  
  AllComb2 <- AllComb
  colnames(AllComb2) <- paste0("nID_fam_plot_G", 1:nGroups)
  
  fill_allcomb <-
    matrix(
      NA,
      nrow = length(Strategies) - (nrow(AllComb2)),
      ncol = nGroups,
      dimnames = list((nrow(AllComb2) + 1):length(Strategies),
                      colnames(AllComb2))
    )
  
  AllComb2 <- rbind(AllComb2, fill_allcomb)
  
  for (i in names(Strategies)) {
    Strategies[[i]][["N_Strategy"]] <- i
  }
  
  # Effective Number (NE) and Genetic Gain with selection (GS)
  GS.NE <- list()
  for (i in seq_along(Strategies)) {
    # GS
    gs <-
      round(mean(Strategies[[i]]$BV, na.rm = T) / mean(Data_Total[[Trait]], na.rm = T) - 1,
            3)
    # NE
    t <- table(as.factor(Strategies[[i]]$Family))
    # Averages of  Trees number for Mother Tree
    n <-
      length(Strategies[[i]][[id]]) / length(levels(as.factor(Strategies[[i]]$Family)))
    # Variance of Trees number for Mother Tree
    VarID <-
      sum((t - (n)) ^ 2) / (length(Strategies[[i]][[id]]) - 1)
    nf <- length(unique(Strategies[[i]]$Family))
    ne <- round((4 * nf * n) / (n + 3 + (VarID / n)), 1)
    gs.ne <-
      data.frame(
        GS = gs,
        `GS(%)` = scales::percent(gs, accuracy = 0.1),
        NE = ne
      ) %>%
      `colnames<-`(., c("GS", "GS(%)", "NE"))
    gs.ne$Trait <- Trait
    gs.ne$N_Strategy <- i
    gs.ne$N_Progeny <- nrow(Strategies[[i]])
    GS.NE[[i]] <-
      gs.ne[, c("N_Strategy", "N_Progeny", "Trait", "GS", "GS(%)", "NE")]
    GS.NE[[i]]$Strategy <- NA
  }
  
  # Joining BW and BI strategies: Pre final table
  for (i in 1:length(G.BW)) {
    GS.NE[[i]]$Strategy <- "BW"
  }
  for (i in (length(G.BW) + 1):(length(Strategies))) {
    GS.NE[[i]]$Strategy <- "BI"
  }
  GS.NE <- bind_rows(GS.NE)
  #GS.NE
  
  # Summarise thinning strategies
  suppressMessages(
    nFam_Strategy <- bind_rows(Strategies) |>
      data.frame() |>
      mutate(N_Strategy = as.numeric(N_Strategy)) |>
      group_by(N_Strategy, Group) |>
      summarise(nFam = n_distinct(Family)) |>
      pivot_wider(
        id_cols = N_Strategy,
        names_from = "Group",
        values_from = nFam,
        names_prefix = "nFam.",
        values_fill = 0
      )
  )
  
  BV_fam <- left_join(BV_fam,
                      Data_Total[match(BV_fam$Family, Data_Total$Family), ] |>
                        dplyr::select(Family, Group),
                      by = "Family")
  
  # Final table with thinning strategies
  Thinning <- left_join(GS.NE, nFam_Strategy, by = "N_Strategy") |>
    cbind(AllComb2)
  
  Thinning <- left_join(
    Thinning[, c("N_Strategy", "N_Progeny")] |>
      mutate(`N_Progeny(%)` = scales::percent(N_Progeny / nrow(Data_Total), accuracy = 0.1)),
    Thinning,
    by = c("N_Strategy", "N_Progeny")
  )
  
  if (save_table_xlsx == TRUE) {
    require(openxlsx)
    write.xlsx(
      Thinning,
      file = paste0(
        "Thinning_Strategies_for_",
        Trait,
        "_trait_",
        Sys.Date(),
        ".xlsx"
      ),
      overwrite = T
    )
    write.xlsx(
      Strategies,
      file = paste0(
        length(Strategies),
        "strategies_for_",
        nGroups,
        "_groups_",
        Trait,
        "_trait_",
        Sys.Date(),
        ".xlsx"
      ),
      overwrite = T
    )
    write.xlsx(
      BV_fam,
      file = paste0(
        "Summary_of_families",
        "_for_",
        nGroups,
        "_groups_",
        Trait,
        "_trait_",
        Sys.Date(),
        ".xlsx"
      ),
      overwrite = T
    )
  }
  
  if (save_plot_rank == TRUE) {
    ggsave(filename = "families_rank.tiff", plot = families_rank)
  }
  
  return(
    list(
      Thinning = Thinning,
      Strategies = Strategies,
      BV_Fam = BV_fam,
      ggPlot_families_rank = families_rank
    )
  )
}
