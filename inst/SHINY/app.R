# app.R
library(shiny)
library(ggplot2)

# --------- Log-likelihoods (parameterized to be well-behaved for optim) ---------
ll_normal <- function(params, x) {
  mu <- params[1]; sigma <- params[2]
  if (sigma <= 0) return(-Inf)
  sum(dnorm(x, mean = mu, sd = sigma, log = TRUE))
}

ll_exponential <- function(params, x) {
  rate <- params[1]
  if (rate <= 0 || any(x < 0)) return(-Inf)
  sum(dexp(x, rate = rate, log = TRUE))
}

ll_poisson <- function(params, x) {
  lambda <- params[1]
  if (lambda <= 0 || any(x < 0) || any(x != floor(x))) return(-Inf)
  sum(dpois(x, lambda = lambda, log = TRUE))
}

ll_bernoulli <- function(params, x) {
  p <- params[1]
  if (p <= 0 || p >= 1 || any(!(x %in% c(0,1)))) return(-Inf)
  sum(dbinom(x, size = 1, prob = p, log = TRUE))
}

ll_gamma <- function(params, x) {
  shape <- params[1]; rate <- params[2]
  if (shape <= 0 || rate <= 0 || any(x <= 0)) return(-Inf)
  sum(dgamma(x, shape = shape, rate = rate, log = TRUE))
}

# --------- Simulation helpers ---------
simulate_data <- function(dist, n, params, seed) {
  set.seed(seed)
  switch(dist,
         "Normal"     = rnorm(n, mean = params$mu, sd = params$sigma),
         "Exponential"= rexp(n, rate = params$rate),
         "Poisson"    = rpois(n, lambda = params$lambda),
         "Bernoulli"  = rbinom(n, size = 1, prob = params$p),
         "Gamma"      = rgamma(n, shape = params$shape, rate = params$rate)
  )
}

# --------- MLE via optim (maximize log-likelihood) ---------
mle_fit <- function(dist, x, init) {
  fn <- switch(dist,
               "Normal"      = function(p) -ll_normal(p, x),
               "Exponential" = function(p) -ll_exponential(p, x),
               "Poisson"     = function(p) -ll_poisson(p, x),
               "Bernoulli"   = function(p) -ll_bernoulli(p, x),
               "Gamma"       = function(p) -ll_gamma(p, x)
  )
  
  # Box constraints to keep parameters valid
  lower <- switch(dist,
                  "Normal"      = c(-Inf, 1e-6),
                  "Exponential" = c(1e-6),
                  "Poisson"     = c(1e-6),
                  "Bernoulli"   = c(1e-6),
                  "Gamma"       = c(1e-6, 1e-6)
  )
  upper <- switch(dist,
                  "Normal"      = c(Inf, Inf),
                  "Exponential" = c(Inf),
                  "Poisson"     = c(Inf),
                  "Bernoulli"   = c(1 - 1e-6),
                  "Gamma"       = c(Inf, Inf)
  )
  
  res <- optim(par = init, fn = fn, method = "L-BFGS-B", lower = lower, upper = upper)
  list(par = res$par, value = res$value, converged = res$convergence == 0)
}

# --------- Likelihood grids for plotting ---------
grid_ll_1d <- function(dist, x, param_name, grid_vals, fixed_params = NULL) {
  ll <- sapply(grid_vals, function(v) {
    params <- switch(dist,
                     "Normal"      = c(if (param_name == "mu") v else fixed_params["mu"],
                                       if (param_name == "sigma") v else fixed_params["sigma"]),
                     "Exponential" = c(v),
                     "Poisson"     = c(v),
                     "Bernoulli"   = c(v),
                     "Gamma"       = c(if (param_name == "shape") v else fixed_params["shape"],
                                       if (param_name == "rate")  v else fixed_params["rate"])
    )
    switch(dist,
           "Normal"      = ll_normal(params, x),
           "Exponential" = ll_exponential(params, x),
           "Poisson"     = ll_poisson(params, x),
           "Bernoulli"   = ll_bernoulli(params, x),
           "Gamma"       = ll_gamma(params, x)
    )
  })
  data.frame(param = grid_vals, loglik = ll)
}

grid_ll_2d <- function(dist, x, p1_name, p1_vals, p2_name, p2_vals) {
  df <- expand.grid(p1 = p1_vals, p2 = p2_vals)
  df$loglik <- apply(df, 1, function(row) {
    params <- switch(dist,
                     "Normal" = c(mu = row["p1"], sigma = row["p2"]),
                     "Gamma"  = c(shape = row["p1"], rate = row["p2"])
    )
    switch(dist,
           "Normal" = ll_normal(params, x),
           "Gamma"  = ll_gamma(params, x)
    )
  })
  names(df) <- c(p1_name, p2_name, "loglik")
  df
}

ui <- fluidPage(
  titlePanel("MLE Explorer: Five univariate distributions"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dist", "Distribution", choices = c("Normal", "Exponential", "Poisson", "Bernoulli", "Gamma")),
      numericInput("n", "Sample size", value = 200, min = 10, max = 5000, step = 10),
      numericInput("seed", "Random seed", value = 123, min = 1, max = 1e6, step = 1),
      hr(),
      uiOutput("param_controls"),
      hr(),
      checkboxInput("show_data", "Show sample summary", TRUE),
      checkboxInput("show_density", "Show data density/histogram", TRUE),
      hr(),
      uiOutput("grid_controls")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Estimates",
                 h4("Parameter estimates (MLE vs true)"),
                 tableOutput("est_table"),
                 verbatimTextOutput("optim_status")
        ),
        tabPanel("Likelihood plot",
                 h4("Log-likelihood visualization"),
                 plotOutput("ll_plot", height = "420px"),
                 helpText("For two-parameter models, a contour is shown. For one-parameter, a curve is shown.")
        ),
        tabPanel("Data",
                 plotOutput("data_plot", height = "420px"),
                 verbatimTextOutput("data_summary")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Dynamic parameter controls based on distribution
  output$param_controls <- renderUI({
    switch(input$dist,
           "Normal" = tagList(
             sliderInput("mu", "True mean (mu)", min = -5, max = 5, value = 0, step = 0.1),
             sliderInput("sigma", "True sd (sigma)", min = 0.2, max = 5, value = 1, step = 0.1)
           ),
           "Exponential" = sliderInput("rate", "True rate", min = 0.1, max = 5, value = 1, step = 0.1),
           "Poisson" = sliderInput("lambda", "True lambda", min = 0.1, max = 10, value = 2, step = 0.1),
           "Bernoulli" = sliderInput("p", "True p", min = 0.05, max = 0.95, value = 0.4, step = 0.05),
           "Gamma" = tagList(
             sliderInput("shape", "True shape", min = 0.5, max = 5, value = 2, step = 0.1),
             sliderInput("rate_g", "True rate", min = 0.2, max = 5, value = 1, step = 0.1)
           )
    )
  })
  
  # Likelihood grid controls
  output$grid_controls <- renderUI({
    switch(input$dist,
           "Normal" = tagList(
             selectInput("normal_view", "Likelihood view", choices = c("Curve: mu" = "mu", "Curve: sigma" = "sigma", "Contour: (mu, sigma)" = "contour")),
             sliderInput("mu_grid", "Grid for mu", min = -5, max = 5, value = c(-2, 2), step = 0.1),
             sliderInput("sigma_grid", "Grid for sigma", min = 0.2, max = 5, value = c(0.5, 2.5), step = 0.1)
           ),
           "Exponential" = sliderInput("rate_grid", "Grid for rate", min = 0.1, max = 5, value = c(0.2, 2.5), step = 0.1),
           "Poisson" = sliderInput("lambda_grid", "Grid for lambda", min = 0.1, max = 10, value = c(0.2, 6), step = 0.1),
           "Bernoulli" = sliderInput("p_grid", "Grid for p", min = 0.01, max = 0.99, value = c(0.05, 0.95), step = 0.01),
           "Gamma" = tagList(
             selectInput("gamma_view", "Likelihood view", choices = c("Curve: shape" = "shape", "Curve: rate" = "rate", "Contour: (shape, rate)" = "contour")),
             sliderInput("shape_grid", "Grid for shape", min = 0.5, max = 5, value = c(0.8, 3.5), step = 0.1),
             sliderInput("rate_grid_g", "Grid for rate", min = 0.2, max = 5, value = c(0.3, 3), step = 0.1)
           )
    )
  })
  
  # Reactive: true parameters
  true_params <- reactive({
    switch(input$dist,
           "Normal"     = list(mu = input$mu, sigma = input$sigma),
           "Exponential"= list(rate = input$rate),
           "Poisson"    = list(lambda = input$lambda),
           "Bernoulli"  = list(p = input$p),
           "Gamma"      = list(shape = input$shape, rate = input$rate_g)
    )
  })
  
  # Reactive: data
  x <- reactive({
    simulate_data(input$dist, input$n, true_params(), input$seed)
  })
  
  # Initial guesses (simple, data-driven)
  init_guess <- reactive({
    xx <- x()
    switch(input$dist,
           "Normal" = c(mean(xx), sd(xx)),
           "Exponential" = c(1 / mean(xx)),
           "Poisson" = c(mean(xx)),
           "Bernoulli" = c(mean(xx)),
           "Gamma" = {
             m <- mean(xx); v <- var(xx)
             # Method of moments as starting point: shape ~ m^2 / v, rate ~ shape / m
             shape0 <- max(1e-3, m^2 / max(v, 1e-6))
             rate0  <- max(1e-3, shape0 / m)
             c(shape0, rate0)
           }
    )
  })
  
  # Fit MLE
  mle <- reactive({
    mle_fit(input$dist, x(), init_guess())
  })
  
  # Estimates table
  output$est_table <- renderTable({
    est <- mle()$par
    tp  <- true_params()
    if (input$dist == "Normal") {
      data.frame(Parameter = c("mu", "sigma"),
                 True = c(tp$mu, tp$sigma),
                 MLE  = round(est, 4))
    } else if (input$dist == "Exponential") {
      data.frame(Parameter = "rate", True = tp$rate, MLE = round(est[1], 4))
    } else if (input$dist == "Poisson") {
      data.frame(Parameter = "lambda", True = tp$lambda, MLE = round(est[1], 4))
    } else if (input$dist == "Bernoulli") {
      data.frame(Parameter = "p", True = tp$p, MLE = round(est[1], 4))
    } else {
      data.frame(Parameter = c("shape", "rate"),
                 True = c(tp$shape, tp$rate),
                 MLE  = round(est, 4))
    }
  }, digits = 4)
  
  output$optim_status <- renderText({
    est <- mle()
    paste0("Converged: ", est$converged, " | Negative log-likelihood at optimum: ", round(est$value, 4))
  })
  
  # Likelihood plots
  output$ll_plot <- renderPlot({
    xx <- x()
    est <- mle()$par
    tp <- true_params()
    
    if (input$dist == "Normal") {
      mu_seq <- seq(input$mu_grid[1], input$mu_grid[2], length.out = 200)
      sigma_seq <- seq(input$sigma_grid[1], input$sigma_grid[2], length.out = 200)
      
      if (input$normal_view == "mu") {
        df <- grid_ll_1d("Normal", xx, "mu", mu_seq, fixed_params = c(mu = est[1], sigma = est[2]))
        g <- ggplot(df, aes(param, loglik)) +
          geom_line(color = "#2C7FB8") +
          geom_vline(xintercept = est[1], color = "#2C7FB8", linetype = "dashed") +
          geom_vline(xintercept = tp$mu, color = "#D95F0E") +
          labs(x = "mu", y = "log-likelihood", title = "Normal: log-likelihood vs mu (sigma fixed at MLE)")
      } else if (input$normal_view == "sigma") {
        df <- grid_ll_1d("Normal", xx, "sigma", sigma_seq, fixed_params = c(mu = est[1], sigma = est[2]))
        g <- ggplot(df, aes(param, loglik)) +
          geom_line(color = "#2C7FB8") +
          geom_vline(xintercept = est[2], color = "#2C7FB8", linetype = "dashed") +
          geom_vline(xintercept = tp$sigma, color = "#D95F0E") +
          labs(x = "sigma", y = "log-likelihood", title = "Normal: log-likelihood vs sigma (mu fixed at MLE)")
      } else {
        df <- grid_ll_2d("Normal", xx, "mu", mu_seq, "sigma", sigma_seq)
        g <- ggplot(df, aes(mu, sigma, z = loglik)) +
          geom_contour_filled(bins = 12) +
          geom_point(aes(x = est[1], y = est[2]), color = "#2C7FB8", size = 3) +
          geom_point(aes(x = tp$mu, y = tp$sigma), color = "#D95F0E", size = 3) +
          labs(x = "mu", y = "sigma", title = "Normal: log-likelihood contour")
      }
      print(g)
    }
    
    if (input$dist == "Exponential") {
      rate_seq <- seq(input$rate_grid[1], input$rate_grid[2], length.out = 200)
      df <- grid_ll_1d("Exponential", xx, "rate", rate_seq)
      ggplot(df, aes(param, loglik)) +
        geom_line(color = "#2C7FB8") +
        geom_vline(xintercept = est[1], color = "#2C7FB8", linetype = "dashed") +
        geom_vline(xintercept = tp$rate, color = "#D95F0E") +
        labs(x = "rate", y = "log-likelihood", title = "Exponential: log-likelihood vs rate")
    }
    
    if (input$dist == "Poisson") {
      lambda_seq <- seq(input$lambda_grid[1], input$lambda_grid[2], length.out = 200)
      df <- grid_ll_1d("Poisson", xx, "lambda", lambda_seq)
      ggplot(df, aes(param, loglik)) +
        geom_line(color = "#2C7FB8") +
        geom_vline(xintercept = est[1], color = "#2C7FB8", linetype = "dashed") +
        geom_vline(xintercept = tp$lambda, color = "#D95F0E") +
        labs(x = "lambda", y = "log-likelihood", title = "Poisson: log-likelihood vs lambda")
    }
    
    if (input$dist == "Bernoulli") {
      p_seq <- seq(input$p_grid[1], input$p_grid[2], length.out = 200)
      df <- grid_ll_1d("Bernoulli", xx, "p", p_seq)
      ggplot(df, aes(param, loglik)) +
        geom_line(color = "#2C7FB8") +
        geom_vline(xintercept = est[1], color = "#2C7FB8", linetype = "dashed") +
        geom_vline(xintercept = tp$p, color = "#D95F0E") +
        labs(x = "p", y = "log-likelihood", title = "Bernoulli: log-likelihood vs p")
    }
    
    if (input$dist == "Gamma") {
      shape_seq <- seq(input$shape_grid[1], input$shape_grid[2], length.out = 150)
      rate_seq  <- seq(input$rate_grid_g[1], input$rate_grid_g[2], length.out = 150)
      
      if (input$gamma_view == "shape") {
        df <- grid_ll_1d("Gamma", xx, "shape", shape_seq, fixed_params = c(shape = est[1], rate = est[2]))
        ggplot(df, aes(param, loglik)) +
          geom_line(color = "#2C7FB8") +
          geom_vline(xintercept = est[1], color = "#2C7FB8", linetype = "dashed") +
          geom_vline(xintercept = tp$shape, color = "#D95F0E") +
          labs(x = "shape", y = "log-likelihood", title = "Gamma: log-likelihood vs shape (rate fixed at MLE)")
      } else if (input$gamma_view == "rate") {
        df <- grid_ll_1d("Gamma", xx, "rate", rate_seq, fixed_params = c(shape = est[1], rate = est[2]))
        ggplot(df, aes(param, loglik)) +
          geom_line(color = "#2C7FB8") +
          geom_vline(xintercept = est[2], color = "#2C7FB8", linetype = "dashed") +
          geom_vline(xintercept = tp$rate, color = "#D95F0E") +
          labs(x = "rate", y = "log-likelihood", title = "Gamma: log-likelihood vs rate (shape fixed at MLE)")
      } else {
        df <- grid_ll_2d("Gamma", xx, "shape", shape_seq, "rate", rate_seq)
        ggplot(df, aes(shape, rate, z = loglik)) +
          geom_contour_filled(bins = 12) +
          geom_point(aes(x = est[1], y = est[2]), color = "#2C7FB8", size = 3) +
          geom_point(aes(x = tp$shape, y = tp$rate), color = "#D95F0E", size = 3) +
          labs(x = "shape", y = "rate", title = "Gamma: log-likelihood contour")
      }
    }
  })
  
  # Data plots and summary
  output$data_plot <- renderPlot({
    xx <- x()
    if (input$show_density) {
      if (input$dist %in% c("Normal", "Exponential", "Gamma")) {
        ggplot(data.frame(x = xx), aes(x)) +
          geom_histogram(aes(y = ..density..), bins = 40, fill = "#9ECAE1", color = "white") +
          geom_density(color = "#2C7FB8", size = 1.2) +
          labs(title = paste(input$dist, "data"), x = "x", y = "density")
      } else {
        ggplot(data.frame(x = xx), aes(x)) +
          geom_histogram(binwidth = 1, fill = "#9ECAE1", color = "white") +
          labs(title = paste(input$dist, "data (counts)"), x = "x", y = "count")
      }
    }
  })
  
  output$data_summary <- renderText({
    if (!input$show_data) return("")
    xx <- x()
    paste(
      "n =", length(xx),
      "| mean =", round(mean(xx), 4),
      "| sd =", round(sd(xx), 4),
      "| min =", round(min(xx), 4),
      "| max =", round(max(xx), 4)
    )
  })
}

shinyApp(ui, server)

