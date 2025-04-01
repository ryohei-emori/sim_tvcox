# 必要なパッケージのロード
#install.packages("msm")
#install.packages("survival")
require(msm)
require(survival)

# CREATING g() AND g^-1()
g.inv <- sqrt
g <- function(x) { x^2 }

# CREATING THE TIME SCALE AND TRANSFORMED TIME SCALE
t <- 0:199
t.diff <- (t[-1] - t[1:(length(t) - 1)])[-(length(t) - 1)]
g.inv.t <- g.inv(t)
g.inv.t.diff <- (g.inv(t[-1]) - g.inv(t[1:(length(t) - 1)]))[-(length(t) - 1)]

# CREATING THE BOUNDS OF TRUNCATION
t.max <- 150
t.min <- 10
g.inv.t.max <- g.inv(t.max)
g.inv.t.min <- g.inv(t.min)

# DATA GENERATING PROCESS FOR COVARIATES
B <- function(N, m, M) { matrix(runif(N * 10, m, M), ncol = 10) }

# BETA COEFFICIENTS
b <- rnorm(10, 0, 1)  # 10個の係数パラメータ

# NUMBER OF OBSERVATIONS
n <- 1000

# CREATING DATA VECTOR
z.list <- list()
for (i in 1:n) {
  z <- B(length(t), -0.5, 0.5)
  z.list[[i]] <- cbind(z, exp(z %*% b))
}

# GENERATING DATA USING ACCEPT-REJECT METHOD
k <- function(x, m, M, rates, t){
  ifelse(x <= m | x >= M, 0, dpexp(x, rates, t))
}

gen.y <- function(x) {
  x1 <- x[, 11]
  d <- ppexp(g.inv.t.max, x1, g.inv.t) - ppexp(g.inv.t.min, x1, g.inv.t)
  M <- 1 / d
  r <- 60
  repeat{
    y <- rpexp(r, x1, g.inv.t)
    u <- runif(r)
    t <- M * ((k(y, g.inv.t.min, g.inv.t.max, x1, g.inv.t) / d / dpexp(y, x1, g.inv.t)))
    y <- y[u <= t][1]
    if (!is.na(y)) break
  }
  y
}

y <- sapply(z.list, gen.y)
g.y <- g(y)

# CREATING CENSORING INDICATOR
prop.cen <- 0.5
d <- sample(0:1, n, replace = TRUE, prob = c(prop.cen, 1 - prop.cen))

# linear dynamic step :発症確率が低い人ほどタイムステップは短い
get_timestep_length <- function(rate_parameter, rate_parameters) {
  # レートパラメータに基づいてタイムステップの長さを決定
  # 例えば、線形スケーリングを使用する場合
  min_timestep <- 1
  max_timestep <- 10
  scaled_length <- min_timestep + (max_timestep - min_timestep) * (rate_parameter / max(rate_parameters))
  return(ceiling(scaled_length))  # 整数に丸める
}

# 逆U字
get_timestep_length <- function(rate_parameter, rate_parameters) {
  median_rate <- median(rate_parameters)
  difference_squared <- (rate_parameter - median_rate)^2
  # タイムステップの長さを計算（中央値で最大、そこから離れるほど短くなる）
  max_timestep <- 10
  min_timestep <- 1
  range <- max_timestep - min_timestep
  # 中央値からの差が大きいほど、タイムステップは短くなる
  timestep_length <- max_timestep - (difference_squared / max(rate_parameters - median_rate)^2) * range
  # タイムステップの長さを整数に丸める
  return(ceiling(timestep_length))
}


# 逆U字＋noise
get_timestep_length <- function(rate_parameter, rate_parameters) {
  median_rate <- median(rate_parameters)
  difference_squared <- (rate_parameter - median_rate)^2
  max_timestep <- 10
  min_timestep <- 1
  range <- max_timestep - min_timestep
  # ノイズを加える (ここでは標準偏差0.5の正規分布を使用)
  noise <- rnorm(1, mean = 0, sd = 0.5)
  # ノイズを加えてタイムステップの長さを計算
  timestep_length <- max_timestep - (difference_squared / max(rate_parameters - median_rate)^2) * range + noise
  # タイムステップの長さを1以上10以下に制限
  return(max(min(ceiling(timestep_length), max_timestep), min_timestep))
}


# CREATING DATASET WITH INDIVIDUAL RATE-DEPENDENT TIMESTEPS
data <- NULL
rate_parameters <- sapply(z.list, function(x) x[1, 11])
weights <- numeric(n) # Initialize the weights vector
for (i in 1:n) {
  rate_param <- z.list[[i]][1, 11]
  timestep_length <- get_timestep_length(rate_param, rate_parameters)
  
  # 個人の観測値の最大値（g.y_val）を取得し、タイムステップのシーケンスを生成
  g.y_val <- ceiling(g(y[i]))
  time_steps <- seq(1, g.y_val, by = timestep_length)
  if (max(time_steps) < g.y_val) {
    time_steps <- c(time_steps, g.y_val) # 最後の観測値を確実に含める
  }
  
  weights <- g.y_val / timestep_length
  
  # 共変量の行を time_steps に応じて選択
  z_temp <- z.list[[i]][, 1:10]  # すべての行と最初の10の共変量
  z_step <- z_temp[time_steps, ]  # time_steps に対応する行だけを選択
  
  # ID、時間、およびイベント指標のデータを作成
  id.temp <- rep(i, length(time_steps))
  time.temp <- time_steps
  time0.temp <- c(0, time_steps[-length(time_steps)])
  d.temp <- c(rep(0, length(time_steps) - 1), d[i])
  
  data.temp <- cbind(id = id.temp, t = time.temp, t0 = time0.temp, d = d.temp, z_step, w=weights)
  data <- rbind(data, data.temp)
}

colnames(data) <- c("id", "t", "t0", "d", paste0("z", 1:10), "w")
data = data.frame(data)

# APPLYING COX MODEL
model <- coxph(Surv(t0, t, d) ~ . -id -w, data = data)
schoenfeld <- cox.zph(model, transform = "identity")

# RESULT
#data
summary(model)
schoenfeld

# APPLYING COX MODEL with Weights
modelw <- coxph(Surv(t0, t, d) ~ . -id -w, data = data, weights = w)
schoenfeldw <- cox.zph(modelw, transform = "identity")

# RESULT
summary(modelw)
schoenfeld

# COMPARING REGULAR VISITS LONGITUDINAL DATA
# CREATING DATASET
data1 <- NULL
for (i in 1:n) {
  id.temp <- rep(i, ceiling(g.y[i]))
  time.temp <- c(1:ceiling(g.y[i]))
  time0.temp <- 0:(ceiling(g.y[i]) - 1)
  d.temp <- c(rep(0, length(time.temp) - 1), d[i])
  z.temp <- z.list[[i]][1:(ceiling(g.y[i])), 1:10]
  data1.temp <- cbind(id.temp, time.temp, time0.temp, d.temp, z.temp)
  data1 <- rbind(data1, data1.temp)
}
colnames(data1) <- c("id", "t", "t0", "d", paste0("z", 1:10))
data1 <- data.frame(data1)

# APPLYING COX MODEL
model1 <- coxph(Surv(t0, t, d) ~ . -id -t -t0, data = data1, id=id)
schoenfeld1 <- cox.zph(model1, transform = "identity")

# RESULT
data1
summary(model1)
schoenfeld1


# COMPARING SINGLE POINT (TIME-INVARIANT COVARIATES) DATA
# data1 と同じカラム構造を持つ空のデータフレームを作成
data2 <- data.frame(matrix(ncol = ncol(data1), nrow = 0))
colnames(data2) <- colnames(data1)

# ユニークな ID ごとにループ
unique_ids <- unique(data1$id)
for (id in unique_ids) {
  id_rows <- data1[data1$id == id, ]
  covariates <- id_rows[1, grepl("^z", names(id_rows))]
  # 最後の行の t と d の値を取得
  last_t <- tail(id_rows$t, n = 1)
  last_d <- tail(id_rows$d, n = 1)
  # 新しい行を作成
  new_row <- data.frame(id = id, t = last_t, t0 = 0, d = last_d, covariates)
  # 新しいデータフレームに行を追加
  data2 <- rbind(data2, new_row)
}

# APPLYING COX MODEL
model2 <- coxph(Surv(t0, t, d) ~ . -id -t -t0, data = data2, id=id)
schoenfeld2 <- cox.zph(model2, transform = "identity")

# RESULT
data2
summary(model2)
schoenfeld2


# COMPARING MI for REGULAR TIMESTEP DATA
# install.packages("mice")
library(mice)
library(dplyr)

# regular_dataの初期化
regular_data <- NULL
for (i in 1:n) {
  max_t <- ceiling(g(y[i]))
  id_data <- expand.grid(id = i, t = 1:max_t)
  id_data <- id_data %>%
    mutate(t0 = lag(t, default = 0), d = 0) %>%
    arrange(id, t)
  
  # 各共変量をNAで初期化
  for(j in 1:10) {
    id_data[paste0("z", j)] <- NA
  }
  regular_data <- rbind(regular_data, id_data)
}

# `data` から各 `id` の最後の `d` の値を取得
last_d_values <- data %>%
  group_by(id) %>%
  summarize(last_d = last(d)) %>%
  ungroup()

# `regular_data` の各 `id` の最後の行に `last_d` の値をマッピング
regular_data <- regular_data %>%
  group_by(id) %>%
  mutate(
    d = ifelse(row_number() == n(), last_d_values$last_d[match(id, last_d_values$id)], 0)
  ) %>%
  ungroup()

# `data` から共変量 `z1` から `z10` の値を `regular_data` にマッピング
regular_data <- merge(regular_data, data %>% select(id, t, z1:z10), by = c("id", "t"), all.x = TRUE)
regular_data <- regular_data %>% select(-matches("\\.x$"))

# 欠損値をmiceで補完する
mice_mod <- mice(regular_data[, !(names(regular_data) %in% c("id", "t", "t0", "d"))], m = 5, method = 'pmm', seed = 500)

# 補完された規則的な間隔のレコードのみを抽出する
completed_data <- complete(mice_mod, action = 1)
regular_completed_data <- cbind(regular_data[, c("id", "t", "t0", "d")], completed_data)

# APPLYING COX MODEL
model3 <- coxph(Surv(t0, t, d) ~ . -id, data = regular_completed_data)
schoenfeld3 <- cox.zph(model3, transform = "identity")

# RESULT
regular_completed_data
summary(model3)
schoenfeld3



# CREATING DATASET WITH RANDOM DYNAMIC TIMESTEPS
data4 <- NULL
get_random_timestep_length <- function() {
  sample(1:10, 1, replace = TRUE)
}
for (i in 1:n) {
  # タイムステップの長さをランダムに決定
  time_steps <- c(1)
  while(max(time_steps) < ceiling(g(y[i]))){
    next_step <- max(time_steps) + get_random_timestep_length()
    time_steps <- c(time_steps, next_step)
  }
  
  # 共変量の行を time_steps に応じて選択
  z_temp <- z.list[[i]][, 1:10]  # すべての行と最初の10の共変量
  z_step <- z_temp[time_steps, ]  # time_steps に対応する行だけを選択
  
  id.temp <- rep(i, length(time_steps))
  time.temp <- time_steps
  time0.temp <- c(0, head(time_steps, -1))
  d.temp <- c(rep(0, length(time_steps) - 1), tail(d[i:length(time_steps)], 1))
  
  # すべてのカラムを結合して新しい行を作成
  data.temp <- cbind(id = id.temp, t = time.temp, t0 = time0.temp, d = d.temp, z_step)
  
  # 新しい行をdata4に追加
  data4 <- rbind(data4, data.temp)
}

colnames(data4) <- c("id", "t", "t0", "d", paste0("z", 1:10))
data4 <- data.frame(data4)

# APPLYING COX MODEL
model4 <- coxph(Surv(t0, t, d) ~ . -id, data = data4)
schoenfeld4 <- cox.zph(model4, transform = "identity")

# RESULT
data4
summary(model4)
schoenfeld4


# CREATING DATASET WITH RANDOM INDIVIDUAL TIMESTEPS
data5 <- NULL
get_random_timestep_length <- function() {
  sample(1:10, 1)
}
random_timestep_lengths <- replicate(n, get_random_timestep_length())
weights <- numeric(n) # Initialize the weights vector
for (i in 1:n) {
  # 各個人IDに対してランダムに割り当てられたタイムステップの長さを取得
  timestep_length <- random_timestep_lengths[i]
  
  # 個人の観測値の最大値を取得
  g.y_val <- ceiling(g(y[i]))
  
  weights <- g.y_val / timestep_length
  
  # タイムステップのシーケンスを生成
  time_steps <- seq(1, g.y_val, by = timestep_length)
  if (max(time_steps) < g.y_val) {
    time_steps <- c(time_steps, g.y_val) # 最後の観測値を確実に含める
  }
  
  # 共変量の行を time_steps に応じて選択
  z_temp <- z.list[[i]][, 1:10]  # すべての行と最初の10の共変量
  z_step <- z_temp[time_steps, ]  # time_steps に対応する行だけを選択
  
  id.temp <- rep(i, length(time_steps))
  time.temp <- time_steps
  time0.temp <- c(0, head(time_steps, -1))
  d.temp <- c(rep(0, length(time_steps) - 1), d[i])
  
  # すべてのカラムを結合して新しい行を作成
  data.temp <- cbind(id = id.temp, t = time.temp, t0 = time0.temp, d = d.temp, z_step, w = weights)
  
  # 新しい行をdata6に追加
  data5 <- rbind(data5, data.temp)
}

colnames(data5) <- c("id", "t", "t0", "d", paste0("z", 1:10), "w")
data5 <- data.frame(data5)

# APPLYING COX MODEL
model5 <- coxph(Surv(t0, t, d) ~ . -id -w, data = data5)
schoenfeld5 <- cox.zph(model5, transform = "identity")

# RESULT
#data5
summary(model5)
schoenfeld5

# APPLYING COX MODEL with Weights
model6 <- coxph(Surv(t0, t, d) ~ . -id -w, data = data6, weights = w)
schoenfeld6 <- cox.zph(model6, transform = "identity")

# RESULT
summary(model6)
schoenfeld6