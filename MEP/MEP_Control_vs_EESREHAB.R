# analysis title
# date

#Clear R's brain
rm(list=ls()) 

library(dplyr)

# データ入力
df <- data.frame(
  group = factor(rep(c("Control", "EES"), each = 4),
                 levels = c("Control", "EES")),
  amplitude = c(1, 8, 19, 42,
                11, 11, 18, 20)
)

df

wilcox.test(amplitude ~ group, data = df, exact = FALSE)

library(ggplot2)

# 各群の要約（mean, SD, SEM）
summary_df <- df %>%
  group_by(group) %>%
  summarise(
    mean = mean(amplitude),
    sd   = sd(amplitude),
    n    = n(),
    sem  = sd / sqrt(n)
  )

summary_df

ggplot(summary_df, aes(x = group, y = mean)) +
  geom_col(width = 0.6) +  # 棒（平均）
  geom_errorbar(aes(ymin = mean - sem,
                    ymax = mean + sem),
                width = 0.2) +  # エラーバー（SEM）
  geom_jitter(data = df,
              aes(x = group, y = amplitude),
              width = 0.05, height = 0,
              size = 2) +        # 各ラットの値
  labs(x = "", y = "MEP amplitude") +
  theme_classic(base_size = 14)

# 必要なパッケージ
library(ggplot2)
library(dplyr)

# データ
df <- data.frame(
  group = factor(rep(c("Control", "EESREHAB"), each = 4),
                 levels = c("Control", "EESREHAB")),
  amplitude = c(1, 8, 19, 42,
                11, 11, 18, 20)
)

# 各群の統計量（mean + SEM）
summary_df <- df %>%
  group_by(group) %>%
  summarise(
    mean = mean(amplitude),
    sd   = sd(amplitude),
    n    = n(),
    sem  = sd / sqrt(n)
  )

# 色（アップロード画像と近い色）
control_color <- "#4D4D4D"   # ピンク系
ees_color <- "#E377C2"       # シアン系

# プロット
ggplot(summary_df, aes(x = group, y = mean, fill = group)) +
  geom_col(width = 0.65, color = "black") +
  geom_errorbar(aes(ymin = mean - sem,
                    ymax = mean + sem),
                width = 0.2, size = 1) +
  geom_jitter(data = df,
              aes(x = group, y = amplitude),
              width = 0.08, size = 2, color = "white", stroke = 1.2) +
  scale_fill_manual(values = c("Control" = control_color,
                               "EESREHAB" = ees_color)) +
  labs(y = "MEP amplitude (mV)", x = "") +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14)
  )

ggplot(summary_df, aes(x = group, y = mean, fill = group)) +
  geom_col(width = 0.65, color = "black") +
  geom_errorbar(aes(ymin = mean - sem,
                    ymax = mean + sem),
                width = 0.2, linewidth = 1) +
  geom_jitter(data = df,
              aes(x = group, y = amplitude, fill=group),
              width = 0.08,
              size = 3,
              shape = 21,
              color = "black",
              stroke = 1.2) +
  scale_fill_manual(values = c("Control" = control_color,
                               "EESREHAB" = ees_color)) +
  labs(y = "MEP amplitude (mV)", x = "") +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14)
  )