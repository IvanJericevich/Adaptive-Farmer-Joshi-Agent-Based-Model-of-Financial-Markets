;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Preliminaries ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
patches-own [
  strategy ;; "Chartist" or "Fundamentalist"
  strategy-int ;; 1 if "Chartist" or -1 if "Fundamentalist"
  omega ;; Order placed at the close of each day
  x-new ;; Position to be taken for the next day
  x-current ;; Today's position
  x-previous ;; Previous position
  c ;; Capital assignment (number of shares as a percentage of market capitalization)
  v-current ;; Today's perceived value for each fundamentalist
  v-new ;; The perceived value for the next day
]

globals [
  number-bull ;; The number of agents taking a positive position
  number-bear ;; The number of agents taking a negative position
  log-price
  log-price-previous
  log-returns
  f
  delta ;; Parameter allowing for recruitment dynamics
]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Supplementary function ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report sign [number] ;; Function to obtain the sign of a number
  ifelse (number > 0) [
    report 1
  ][
    ifelse (number < 0) [
      report (-1)
    ][
      report 0
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Initialisation ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to setup
  __clear-all-and-reset-ticks
  ;; Initialise strategies
  let N-chartists ceiling (count patches * chartist-fundamentalist-ratio);floor (count patches / 2)
  let N-fundamentalists count patches - N-chartists
  ask n-of N-chartists patches [
    set strategy "Chartist"
    set strategy-int 1
  ]
  ask patches with [strategy != "Chartist"] [
    set strategy "Fundamentalist"
    set strategy-int (-1)
  ]
  initialisation
  set f 0.00001
  set delta 0
end

to initialisation ;; Initialise price, return and positions by sequentially applying a random walk to the initial price and updating positions (10 days)
  set log-price 50 ;; Initial price
  ask patches [
    ;; Initialise positions and value perceptions
    set c random-float 1
    set x-previous c
    set x-current c * -1
    set x-new c
    set v-current random-normal (mu_eta + log-price) sigma_eta
    set v-new random-normal (mu_eta + log-price) sigma_eta
  ]
  ;; Simulate agents observing 10 days of prices before taking actual positions
  repeat 10 [
    set log-price-previous log-price
    set log-price log-price + random-normal 0 sigma_zeta
    ask patches [
      update-value-perception
      update-positions
    ]
  ]
end

to update-positions ;; Update the positions of all agents for the next day after orders have been submitted
  set x-previous x-current
  set x-current x-new
  let sign-current sign x-current
  set c random-float 1 ;abs (random-normal 1 10) ;; Random capital assignment as a percentage of market capitalization
  let sign-new sign (ifelse-value (strategy = "Chartist") [log-price - log-price-previous] [v-current - log-price]) ;; Agents only care about the direction of price deviations
  ;set x-new ifelse-value (sign-new = sign-current) [x-current] [c * sign-new]
  set x-new x-current + c * sign-new ;; Calculate the position for the next day
end

to update-value-perception ;; Update fundamentalists's perceived value of the stock
  set v-current v-new
  set v-new random-normal (log-price + mu_eta) sigma_eta
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Strategies ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to trend-followers ;; Trend follower submits order, observes history of prices & information and then sets their positions for the next day accordingly
  ask patches [
    if strategy = "Chartist" [
      ;; Calculate today's order
      set omega x-current - x-previous
      update-positions
      ;; Visualisation
      ifelse (set-patches-to = "Positions") [
        if (omega > 0) [set pcolor scale-color blue (omega) 1 0];[set pcolor blue]
        if (omega < 0) [set pcolor scale-color red (omega) -1 0];[set pcolor red]
        if (omega = 0) [set pcolor black]
      ][
        set pcolor orange
      ]
    ]
  ]
end

to value-investors ;; Value investor submits order, observes history of prices & information and then set their positions for the next day accordingly
  ask patches [
    if strategy = "Fundamentalist" [
      ;; Calculate today's order
      set omega x-current - x-previous
      update-positions
      update-value-perception
      ;; Visualisation
      ifelse (set-patches-to = "Positions") [
        if (omega > 0) [set pcolor blue]
        if (omega < 0) [set pcolor red]
        if (omega = 0) [set pcolor black]
      ][
        set pcolor green
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Market clearing mechanism ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to market-maker ;; Market maker bases price formation on the net orders
  ;; Update previous price
  set log-price-previous log-price
  ;; Noise term
  let zeta random-normal 0 sigma_zeta
  ;; Calculate price
  set log-price log-price + (sum [omega] of patches) / lambda + zeta
  ;; Update return
  set log-returns log-price - log-price-previous
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Agent clustering mechanisms ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to herd
  set f (f + (random-normal 0.5 1))
  let A (abs (sin (f)) * 0.03)
  ask patches [
    ifelse ((A * sum [strategy-int] of neighbors4 + random-normal 0 0.05) > 0) [
      set strategy "Chartist"
      set strategy-int 1
    ][
      set strategy "Fundamentalist"
      set strategy-int -1
    ]
  ]
  set chartist-fundamentalist-ratio (count patches with [strategy-int = 1]) / (count patches) ;; round this = round
end

to recruit
  let N count patches
  let some random N
  let epsilon 0.5
  set delta abs sin (delta + 0.01)
  let theta count patches with [strategy = "Chartist"]; + random 5
  let pchartist (1 - theta / N) * (epsilon + (1 - delta) * (theta / (N - 1)))
  let pfundamentalist (theta / N) * (epsilon + (1 - delta) * (N - theta) / (N - 1))
  let uniform random-float 1
  ask one-of patches [
    ask patches in-radius 3 [
      ifelse ((strategy = "Fundamentalist") and (uniform <= pchartist)) [
        set strategy "Chartist"
        set strategy-int 1
      ][
        if ((strategy = "Chartist") and (uniform > pchartist) and (uniform <= (pchartist + pfundamentalist))) [
          set strategy "Fundamentalist"
          set strategy-int (-1)
        ]
      ]
    ]
  ]
  set chartist-fundamentalist-ratio (count patches with [strategy-int = 1]) / (count patches) ;; round this = round
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Implementation ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to go
  ;; Implement strategies
  trend-followers
  value-investors
  ;; Orders are submitted to the market maker
  market-maker
  ;; Clustering of agents
  if (adaptation = "herding") [herd]
  if (adaptation = "recruitment") [recruit]
  ;; Visualisation
  do-plot
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Visualisation ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to do-plot
  ;; Visualise prices and returns
  set-current-plot "Log-Price"
  set-current-plot-pen "Log-Price"
  plot log-price
  set-current-plot "Log-Returns"
  set-current-plot-pen "Log-Return"
  plot log-returns
  ;; Visualise positions
  set number-bull count patches with [omega > 0]
  set number-bear count patches with [omega < 0]
  set-current-plot "Positions"
  set-current-plot-pen "Bulls"
  plot number-bull
  set-current-plot-pen "Bears"
  plot number-bear
  ;; Visualise strategies
  set-current-plot "Strategies"
  set-current-plot-pen "Fundamentalists"
  plot count patches with [strategy-int = -1]
  set-current-plot-pen "Chartists"
  plot count patches with [strategy-int = 1]
end
@#$#@#$#@
GRAPHICS-WINDOW
216
10
554
349
-1
-1
10.0
1
10
1
1
1
0
1
1
1
-16
16
-16
16
0
0
1
ticks
30.0

BUTTON
8
10
98
43
Go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
118
10
207
43
Setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
9
77
208
110
lambda
lambda
0
15
7.5
0.01
1
NIL
HORIZONTAL

SLIDER
9
134
209
167
sigma_zeta
sigma_zeta
0
10
2.5
0.01
1
NIL
HORIZONTAL

PLOT
561
183
863
349
Log-Returns
Time
Log-Return
0.0
10.0
-10.0
10.0
true
false
"" ""
PENS
"Log-Return" 1.0 0 -16777216 true "" ""

PLOT
560
10
863
178
Log-Price
Time
Log-Price
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"Log-Price" 1.0 0 -16777216 true "" ""

PLOT
869
10
1215
178
Positions
Time
Counts
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Bears" 1.0 0 -5298144 true "" ""
"Bulls" 1.0 0 -13345367 true "" ""

SLIDER
9
207
209
240
mu_eta
mu_eta
-50
50
0.0
0.01
1
NIL
HORIZONTAL

SLIDER
9
248
209
281
sigma_eta
sigma_eta
0
10
5.0
0.01
1
NIL
HORIZONTAL

SLIDER
9
318
208
351
chartist-fundamentalist-ratio
chartist-fundamentalist-ratio
0
1
0.5
0.001
1
NIL
HORIZONTAL

PLOT
869
183
1260
349
Strategies
Time
Counts
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Fundamentalists" 1.0 0 -13840069 true "" ""
"Chartists" 1.0 0 -955883 true "" ""

CHOOSER
10
379
207
424
set-patches-to
set-patches-to
"Positions" "Strategies"
0

TEXTBOX
65
56
158
74
Liquidity parameter
11
0.0
1

TEXTBOX
74
115
143
133
Noise process
11
0.0
1

TEXTBOX
67
172
161
200
Fundamentalists' sentiment settings
11
0.0
1

TEXTBOX
65
286
167
314
Ratio of chartists to fundamentalists
11
0.0
1

TEXTBOX
62
357
151
375
Patch visualisation
11
0.0
1

TEXTBOX
227
354
377
372
Agent adaptation technique
11
0.0
1

CHOOSER
240
388
378
433
adaptation
adaptation
"<none>" "herding" "recruitment"
0

@#$#@#$#@
# WHAT IS IT?

This model is an inter-day ABM that uses a market maker based method of price
formation to study the price dynamics induced by two commonly used financial trading strategies - trend following (where agents are called chartists) and value investing (where agents are called fundamentalists). This is an interesting model for understanding daily trading decisions made from closing auction to closing auction in equity markets, as it attempts to model financial market behaviour without the inclusion of agent adaptation. As an extension, however, I consider the case where agents are allowed to switch strategies probabilistically using either a herding or recruiment (following Kirman and Teyssiere (2002)) approach. In the herding approach strategies are favoured according to the periodic degree of importance placed on the positions of each agent's neighbours. On the other hand, in the trader-agent recruitment approach, groups of agents in a neighbourhood come together and try to recruit the minority to the majority strategy in that neighbourhood. The aim here is to assign realistic rules to trading agents such that important fetaures of financial market time-series (known as stylized facts - eg. leptokurtosis and significant auto-correlations of absolute log-returns) can be replicated/recovered.

The market maker determines the price using the market impact function which relates the net of all orders at a given time to prices. Trend followers invest based on the belief that price changes have inertia while value investors believe their perceived value may not be fully reflected in the current price, and that the price will move towards their perceived value (their decisions are based on price deviations from a perceived/subjective fundamental asset value). Each strategy induces price dynamics that characterize its signal processing properties.

# HOW IT WORKS

This model only considers the case where agents trade in a single stock and submit market orders (the model studies only market orders) as a percentage of the stocks market capitalization (U(0, 1)). Individual trader agents are given by the patches who can buy or sell a specified amount of the asset at each time. The number of agents remain fixed throughout the simulation. The model assumes chartists care only about the direction of the change from the lag-1 price to the current price and not about the magnitude. Similarly, fundamentalists care only about the direction of the change from the current price to the perceived value. During each of T discrete simulation days, the following occurs:
1. Each trader agent observes the most recent prices and the information and
submits orders to buy or sell some quantity of an asset to a risk-neutral market maker at day-end.
2. The market maker then fills the requested orders at a newly determined market price based on a closed-form equation aggregating trader agent demands (market-impact function).

More specifically, the simulation begins by initialising prices and positions of all agents for 10 days. Therefater the actual simulation begins - trend followers enter a positive position if the price increased from yesterday and vise versa while value investors make a subject assessment of value with uncertainty (N(mu_eta, sigma_eta)) and then take a positive position if the stock price is below this value (and vice versa). The market maker then sets the next day's price according to the market impact function, afterwhich agents revise their strategy based on the clustering mechanism. The herding clustering mechanism works by sequentially updating the degree of importance an agent gives to his neighbours. This behaviour is periodic - that is, if the importance is high then we observe herding behaviour if it is small then there is a disorganized state.

# HOW TO USE IT

## Agent settings
The values percieved or sentiments adopted by all fundamentalists are determined by a random normal distribution where their average sentiment is controlled by the mu_eta slider (large mu_eta implies positive sentiment while small mu_eta implies negative sentiment) and the uncertainty in value perceptions is controlled by the sigma_eta slider. Additionally, the ratio of chartists to fundamentalists (when agents don't herd) can be adjusted as well. More importantly, the method through which agents cluster is chosen by using the adaptation selector.

## Price settings
The lambda slider represents the liquidity of the market maker and determines the effect that the size of orders have on price formation - that is, the greater the liquidity in the market, the larger the sizes of the orders placed and hence their effect on price. Furthermore, the sigma_zeta slider controls the amount of noise in the price formation process and can be interpreted as corresponding to noise traders who submit orders at random or as random information that affects the market maker’s price setting decisions. Lastly, the initial price may be adjusted as well.

# THINGS TO NOTICE

## Without herding
The most realistic results without herding are obtained when there are almost an equal number of chartists and fundamentalists. A large proportion of either strategy will result in an almost linear price series where with a consistent downward/upward trend. Overall, deviating from a chartist-to-fundamentalist ratio of 0.5 produces very unrealistic results. Notice that, with an equal number of chartists and fundamentalists, when setting mu_eta < 0 (bearish fundamentalist sentiment) then there is a general downward trend in the price series. Similarly, with mu_eta > 0 (bullish fundamentalist sentiment) there is a general upward trend in the price series.

## With herding/recruitment
The most realistic results with herding are obtained when mu_eta = 0 (neutral fundamentalists) and sigma_zeta > 0. Setting mu_eta > 0 will result in an upward trend while mu_eta < 0 will result in a downward trend. Observe the periods where a single strategy dominates and see how the price series is affected.

With herding/recruitment activated (and switching between patch visualisations), notice how the market goes from disorganized states with great diversity of positions to organized behaviour where there are clusters of the same positions regarding the buying or selling of the stock. This then characterises the crashes and booms in the market.

# THINGS TO TRY

## Without herding
1. Set the ratio of chartists to fundamentalists to 1 and then 0 (corresponding to only chartists amd only fundamentalists) and notice the strong linear trends. This is unrealistic but corresponds to the rules assigned to agents.

2. Set the chartist-to-fundamentalist ratio to 0.5 and notice the more realistic price series. Having an equal number of chartists and fundamentalists while setting fundamentalist sentiment to be positive causes a general upward trend. Similarly, a negative sentiment results in a downward trend.

3. Experiment with different degrees of noise.

## With herding/recruitment
1. Using the default parameters, switch between to two patch visualisation selectors to look at the position and strategy dynamics respectively.

2. Try setting fundamentalist sentiment to -25 (bearish - price is above the fundamental price), 0 (neutral - price is at the fundamental price) and 25 (bullish - price is below the fundamental price).

3. Try setting the noise process to zero with equal numbers of chartists and fundamentalists and mu_eta = 0 (neutral fundamentalist sentiment). Although this produces a relatively unrealistic price series, the effect of strategy herding is much more easily seen.

# EXTENDING THE MODEL

* There are many other common strategies that are observed in financial markets which may be considered (eg. contrarians).
* In this cases agents cluster in a very simplified manner. It may be useful for agents to switch their strategies according trading strategy profitability such that many agents adopt the most profitable strategy at each time. Alternatively, their are many other methods for clustering: minority-game, predator-prey, etc.
* This model has not considered the case where agents hold profitable orders/positions or remain neutral for a specified period of time (long/short term). Instead, agents place new orders every day. This can be an unrealistic assumption since placing daily orders would involve large transaction costs.
* It may be useful to allow for agents to die-out when profitability falls below a certain threshold and for new agents to enter. In this way, the number of agents in the system is constantly changing.
* Slighlty more challenging is to apply this model to the case where both market and limit orders are considered.

# RELATED MODELS

Farmer, J. D. and Joshi, S. (2002). The Price Dynamics of Common Trading Strategies. Journal of Economic Behavior & Organization, 49 (2), 149–171.

Alan Kirman and Gilles Teyssiere. Microeconomic models for long memory in the volatility of financial time series. Studies in Nonlinear Dynamics & Econometrics, 5(4), 2002.
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="Standard" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <metric>log-price</metric>
    <metric>log-returns</metric>
    <enumeratedValueSet variable="lambda">
      <value value="7.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-patches-to">
      <value value="&quot;Positions&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="adaptation">
      <value value="&quot;&lt;none&gt;&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chartist-fundamentalist-ratio">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sigma_zeta">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mu_eta">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sigma_eta">
      <value value="5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Herding" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <metric>log-price</metric>
    <metric>log-returns</metric>
    <enumeratedValueSet variable="lambda">
      <value value="7.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-patches-to">
      <value value="&quot;Positions&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="adaptation">
      <value value="&quot;herding&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chartist-fundamentalist-ratio">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sigma_zeta">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mu_eta">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sigma_eta">
      <value value="5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Recruitment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <metric>log-price</metric>
    <metric>log-returns</metric>
    <enumeratedValueSet variable="lambda">
      <value value="7.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-patches-to">
      <value value="&quot;Positions&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chartist-fundamentalist-ratio">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sigma_zeta">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mu_eta">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="adaptation">
      <value value="&quot;recruitment&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sigma_eta">
      <value value="5"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
