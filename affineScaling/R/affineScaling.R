#' Title
#'
#' @param A The constraint Matrix
#' @param x initiadata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAbElEQVR4Xs2RQQrAMAgEfZgf7W9LAguybljJpR3wEse5JOL3ZObDb4x1loDhHbBOFU6i2Ddnw2KNiXcdAXygJlwE8OFVBHDgKrLgSInN4WMe9iXiqIVsTMjH7z/GhNTEibOxQswcYIWYOR/zAjBJfiXh3jZ6AAAAAElFTkSuQmCCl values
#' @param C Objective Function Matrix
#'
#' @return Optimal solutions and objective Function value
#' @export
#'
#' @examples
#' affineScaling(A,x,C)
affineScaling <- function(A,x,C){
  A <- standard(A)
  max_iter <- 100
  X <- diag(x)
  k <- 1
  while(k < max_iter){
    cat("\n________________________________________________ itération", k, "_________________________________________________________\n\n")
    valeurs_approchees(A,X,C)
    couts_reduit(A,X,C)
    direction_translation(A,X,C)
    if(optimal(A,X,C,tolerance = 1e-8) == TRUE){
      cat("\nSolution optimale trouvée à la", k,"e itération\n\n")
      cat("Les solutions sont :", x, "Et la fonction objective est :", sum(t(C)%*%X), "\n\n")
      break
    }
    if(finie(A,X,C) == TRUE){
      cat("\nSolution non finie\n\n")
      break
    }
    calcul_pas(A,X,C)
    #Mise à jour de la valeur de x
    #Affichage
    Valeurs_x <- x
    Couts_reduit <- couts_reduit(A,X,C)
    Direction_translation <- direction_translation(A,X,C)
    Affiche_valeurs <- cbind(Valeurs_x, Couts_reduit, Direction_translation)
    colnames(Affiche_valeurs) <- c("valeurs de x", "couts_reduits", "Direction de translation")
    cat("\n valeurs approchées :", valeurs_approchees(A,X,C), "\n\n")
    print(Affiche_valeurs)
    cat("\n La fonction objective est :", sum(t(C)%*%X), "\n")
    x <- x + calcul_pas(A,X,C)*X%*%direction_translation(A,X,C)
    k <- k+1
    X <- diag(c(x))
  }
  return("Terminé")
}

#Forme standard
standard <- function(M){
  nombre_colonnes <- ncol(M) #nombre de variables
  nombre_lignes <- nrow(M) #nombre de contraintes
  Forme_standard <- cbind(M, diag(nombre_lignes))
  return(Forme_standard)
}

#1ere fonction (valeurs approchées)
valeurs_approchees <- function(A,X,C){
  carree_X <- X%*%X
  produit_1 <- solve((A%*%carree_X%*%t(A)))
  produit_2 <- (A%*%carree_X%*%C)
  produit <- produit_1%*%produit_2
  return(produit)
}

#2e fonction (couts reduits)
couts_reduit <- function(A,X,C){
  produit_1 <- t(A)%*%valeurs_approchees(A,X,C)
  resultat <- C - produit_1
  return(resultat)
}

#3e fonction (direction de translation)
direction_translation <- function(A,X,C){
  resultat <- -(X%*%couts_reduit(A,X,C))
  return(resultat)
}

#4e fonction (vérification de l'optimatilité)
optimal <- function(A,X,C,tolerance=1e-8){
  negatif <- 0
  for (i in couts_reduit(A,X,C)) {
    if (i < 0)
      negatif <-negatif+1
  }
  if(negatif == 0 & sum(X%*%couts_reduit(A,X,C))<tolerance)
    return(TRUE)
  else
    return(FALSE)
}

#5e fonction (Test de finitude)
finie <- function(A,X,C){
  negatif <- 0
  for (i in direction_translation(A,X,C)) {
    if (i < 0)
      negatif <-negatif+1
  }
  if(negatif > 0)
    return(FALSE)
  else
    return(TRUE)
}

#6e fonction (calcul du pas)
calcul_pas <- function(A,X,C){
  alpha <- 0.99
  distance <- 0
  for(i in direction_translation(A,X,C)){
    if(i < 0){
      distance <- c(distance, i)
    }
  }
  min_distance <- min(distance)
  pas <- alpha/-min_distance
  return(pas)
}
