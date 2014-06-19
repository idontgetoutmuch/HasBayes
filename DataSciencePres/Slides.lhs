\documentclass{beamer}
%include polycode.fmt
%options ghci -fglasgow-exts

%format rho = "\rho"
%format y1 = "y_1"
%format y2 = "y_2"
%format var = "(1 - \rho^2)"
%format oldTheta2 = "{\theta}_2"
%format newTheta1 = "\tilde{\theta}_1"
%format newTheta2 = "\tilde{\theta}_2"
%format oldTau    = "{\tau}"
%format newTau    = "\tilde{\tau}"
%format initTau   = "{\tau}_0"
%format newMu     = "\tilde{\mu}"
%format xBar      = "\bar{x}"
%format x2Sum     = "\sum x_i^2"


\usepackage[latin1]{inputenc}
\usepackage{listings}
\usepackage{color}
\usepackage{ulem}

\newcommand {\framedgraphic}[2] {
    \begin{frame}{#1}
        \begin{center}
            \includegraphics[width=\textwidth,height=0.8\textheight,keepaspectratio]{#2}
        \end{center}
    \end{frame}
}

\usetheme{Warsaw}
\setbeamertemplate{footline}{\insertframenumber/\inserttotalframenumber}

\title[Gibbs Sampling]{Gibbs Sampling \\ Some Theory and Some Practice}
\author{Dominic Steinitz}
\date{25th June 2014}
\begin{document}

\begin{frame}
\titlepage
\end{frame}

\section{Introduction}

\subsection{What This Talk is About}

\begin{frame}{What This Talk Is About}

  \begin{itemize}
  \item Introduction to Markov chain methods.
  \item Two Gibbs samplers: theory.
  \item .
  \item .
  \end{itemize}

\end{frame}

\section{Introduction to Markov Chain Methods}

\subsection{The Grasshpper's Problems}

\framedgraphic{Visit Each Hour Equally}{diagrams/ClockEqual12.png}

\begin{frame}[fragile]{Suppose Achieve Goal}

  \uncover<1->{
    $$
    \pi'(n) = p(n \, \vert\, n - 1)\pi(n - 1) + p(n \, \vert\, n + 1)\pi(n + 1)
    $$
  }
  \uncover<2->{
    $$
    \pi'(n) = \frac{1}{2}\frac{1}{N} + \frac{1}{2}\frac{1}{N} = \frac{1}{N}
    $$
  }

\end{frame}

\begin{frame}[fragile]{Testing Convergence}
\begin{verbatim}
ghci> startOnOne
  (1><5)
   [ 1.0, 0.0, 0.0, 0.0, 0.0 ]

ghci> eqProbsMat
  (5><5)
   [ 0.0, 0.5, 0.0, 0.0, 0.5
   , 0.5, 0.0, 0.5, 0.0, 0.0
   , 0.0, 0.5, 0.0, 0.5, 0.0
   , 0.0, 0.0, 0.5, 0.0, 0.5
   , 0.5, 0.0, 0.0, 0.5, 0.0 ]
\end{verbatim}
\end{frame}

\begin{frame}[fragile]{Testing Convergence}
\begin{verbatim}

ghci> take 1 $
      drop 1000 $
      iterate (<> eqProbsMat) startOnOne
  [(1><5)
   [ 0.2, 0.2, 0.2, 0.2, 0.2 ]]
\end{verbatim}
\end{frame}

\framedgraphic{Visit Each Hour Proportionally}{diagrams/ClockProp12.png}

\begin{frame}[fragile]{Suppose Achieve Goal}
  \uncover<1->{
    $$
    \pi'(n) = p(n \, \vert\, n - 1)\pi(n - 1) + p(n \, \vert\, n)\pi(n) + p(n \, \vert\, n + 1)\pi(n + 1)
    $$
  }
  \uncover<2->{
    $$
    \begin{aligned}
      \pi'(4) &= \frac{1}{2}\pi(3) +
      \frac{1}{2}\frac{1}{4}\pi(4) +
      \frac{1}{2}\frac{4}{5}\pi(5) \\
      &= \frac{1}{2}\bigg(\frac{3}{N} + \frac{1}{4}\frac{4}{N} + \frac{4}{5}\frac{5}{N}\bigg) \\
      &= \frac{1}{N}\frac{8}{2} \\
      &= \frac{4}{N} \\
      &= \pi(4)
    \end{aligned}
    $$
  }
\end{frame}

\begin{frame}[fragile]{Testing Convergence}
\begin{verbatim}
    ghci> incProbsMat
    (5><5)
    [ 0.00, 0.500, 0.000, 0.000, 0.5
    , 0.25, 0.250, 0.500, 0.000, 0.0
    , 0.00, 0.333, 0.167, 0.500, 0.0
    , 0.00, 0.000, 0.375, 0.125, 0.5
    , 0.10, 0.000, 0.000, 0.400, 0.5 ]

    ghci> take 1 $
          drop 1000 $
          iterate (<> incProbsMat) startOnOne
    [(1><5)
      [ 6.67e-2, 0.133, 0.199, 0.267, 0.333 ]]
\end{verbatim}
\end{frame}

\begin{frame}[fragile]{Very Simple Markov Chain}
$$
q(i, j) =
\begin{cases}
\mathbb{P} (X_{n+1} = j \,\vert\, X_n = i) = \frac{1}{2} & \text{if } j = i + 1 \mod N \\
\mathbb{P} (X_{n+1} = j \,\vert\, X_n = i) = \frac{1}{2} & \text{if } j = i - 1 \mod N \\
\mathbb{P} (X_{n+1} = j \,\vert\, X_n = i) = 0 & \text{otherwise }
\end{cases}
$$
\end{frame}

\begin{frame}[fragile]{Simple Markov Chain}
  \uncover<1->{
    The grasshopper knows that $\pi(i) = i/N$ so it can calculate
    $\pi(j)/\pi(i) = j/i$ without knowing $N$.
  }
\uncover<2->{
  $$
  p(i, j) =
  \begin{cases}
    q(i,j)\bigg[\frac{\pi(j) q(j,i)}{\pi(i) q(i,j)} \land 1 \bigg] & \text{if } j \ne i \\
    1 - \sum_{k : k \ne i} q(i,k) \bigg[\frac{\pi(k) q(k,i)}{\pi(i) q(i,k)} \land 1 \bigg] & \text{if } j = i
  \end{cases}
  $$
}
\end{frame}

\begin{frame}[fragile]{Simple Markov Chain}
$$
q(i, j) =
\begin{cases}
\frac{1}{2} (\frac{j}{i} \land 1) & \text{if } j \text{ is 1 step clockwise} \\
\frac{1}{2} (\frac{j}{i} \land 1) & \text{if } j \text{ is 1 step anti-clockwise} \\
1 - \frac{1}{2}(\frac{j^c}{i} \land 1) - \frac{1}{2}(\frac{j^a}{i} \land 1) & j = i \& j^a, j^c \text{step (anti-)clockwise} \\
0 & \text{otherwise}
\end{cases}
$$
\end{frame}

\section{Metropolis}

\subsection{Ergodic Theorem}

\begin{frame}{Orientation}
  \begin{itemize}
  \item
    Markov chain theory normally interested in when chain has stationary distribution.
    \pause
  \item
    Here we have a distribution and want to create a Markov chain which has it as stationary distribution.
  \end{itemize}
\end{frame}

\subsection{Metropolis-Hastings Differentiation}

\begin{frame}[fragile]{Ergodic Theorem}
An irreducible, aperiodic and positive recurrent Markov chain has a
unique stationary distribution.

Roughly speaking

\begin{itemize}
  \item
    Irreducible means it is possible to get from any state to any other state.
    \pause
  \item
    Aperiodic means that returning to a state having started at that
    state occurs at irregular times.
    \pause
    \item
      Positive recurrent means that the first time to hit a state is
      finite (for every state and more pedantically except on
      sets of null measure).
\end{itemize}
\end{frame}

\begin{frame}[fragile]{Metropolis-Hastings}
Let $\pi$ be a probability distribution on the state space
$\Omega$ with $\pi(i) > 0$ for all $i$ and let $(Q, \pi_0)$ be an
ergodic Markov chain on $\Omega$ with transition probabilities $q(i,j)> 0$
(the latter condition is slightly stronger than it need be but
we will not need fully general conditions).
\end{frame}

\begin{frame}[fragile]{Metropolis-Hastings}
  Create a new (ergodic) Markov chain with transition probabilities

  $$
  p_{ij} =
  \begin{cases}
    q(i,j)\bigg[\frac{\pi(j) q(j,i)}{\pi(i) q(i,j)} \land 1 \bigg] & \text{if } j \ne i \\
    1 - \sum_{k : k \ne i} q(i,k) \bigg[\frac{\pi(j) q(j,i)}{\pi(i) q(i,j)} \land 1 \bigg] & \text{if } j = i
  \end{cases}
  $$

  where $\land$ takes the maximum of its arguments.
\end{frame}

\begin{frame}[fragile]{Metropolis-Hastings is Ergodic}
$$
\begin{aligned}
\pi(i) q(i,j)\bigg[\frac{\pi(j) q(j, i)}{\pi(i)q(i,j)} \land 1\bigg]
&= \pi(i)q(i,j) \land \pi(j)q(j,i) \\
&= \pi(j)q(j,i)\bigg[\frac{\pi(i) q(i, j)}{\pi(j)q(j,i)} \land 1\bigg]
\end{aligned}
$$

Secondly since we have specified that $(Q, \pi_0)$ is ergodic then
clearly $(P, \pi_0)$ is also ergodic (all the transition probabilities
are $> 0$).

So we know the algorithm will converge to the unique distribution we
specified to provide estimates of values of interest.
\end{frame}

\section{Gibbs Sampling}

\subsection{Random Scan}

\begin{frame}[fragile]{Random Scan}
  Start with $\pi(\theta_1, \theta_2)$.
  $$
  \begin{cases}
    \text{Sample } \theta_1^{(i+1)} \sim \pi(\theta_1 \,\big\vert\, \theta_2^{(i)}) & \text{with probability } \frac{1}{2} \\
    \text{Sample } \theta_2^{(i+1)} \sim \pi(\theta_2 \,\big\vert\, \theta_1^{(i)}) & \text{with probability } \frac{1}{2}
  \end{cases}
  $$

  The transition density kernel is then given by

  $$
  \begin{aligned}
  q\big(\boldsymbol{\theta}^{(i+1)}, \boldsymbol{\theta}^{(i)}\big) &=
  \frac{1}{2}\pi(\theta_1^{(i+1)} \,\big\vert\, \theta_2^{(i)})\delta({\theta_2^{(i)},\theta_2^{(i+1)}}) \\
  &+
  \frac{1}{2}\pi(\theta_2^{(i+1)} \,\big\vert\, \theta_1^{(i)})\delta({\theta_1^{(i)},\theta_1^{(i+1)}})
  \end{aligned}
  $$

  where $\delta$ is the Dirac delta function.
\end{frame}

\begin{frame}[fragile]{Detailed Balance}
  $$
  \begin{aligned}
    & \pi(\theta_1, \theta_2)
    \bigg[
      \frac{1}{2}\pi(\theta_1' \,\big\vert\, \theta_2)\delta({\theta_2,\theta_2'}) +
      \frac{1}{2}\pi(\theta_2' \,\big\vert\, \theta_1)\delta({\theta_1,\theta_1'})\bigg] = \\
    & \frac{1}{2}\bigg[\pi(\theta_1, \theta_2)
      \pi(\theta_1' \,\big\vert\, \theta_2)\delta({\theta_2,\theta_2'}) +
      \ldots
      \bigg] = \\
    & \frac{1}{2}\bigg[\pi(\theta_1, \theta_2')
      \pi(\theta_1' \,\big\vert\, \theta_2)\delta({\theta_2,\theta_2'}) +
      \ldots
      \bigg] = \\
    & \frac{1}{2}\bigg[
      \pi(\theta_2')\pi(\theta_1 \,\big\vert\, \theta_2')
      \pi(\theta_1' \,\big\vert\, \theta_2)\delta({\theta_2,\theta_2'}) +
      \ldots
      \bigg] = \\
    & \frac{1}{2}\bigg[
      \pi(\theta_1', \theta_2')\pi(\theta_1 \,\big\vert\, \theta_2')
      \delta({\theta_2',\theta_2}) +
      \ldots
      \bigg] = \\
    & \pi(\theta_1', \theta_2')\bigg[
      \frac{1}{2}\pi(\theta_1 \,\big\vert\, \theta_2')
      \delta({\theta_2',\theta_2}) +
      \frac{1}{2}\pi(\theta_2 \,\big\vert\, \theta_1')
      \delta({\theta_1',\theta_1})
      \bigg]
  \end{aligned}
  $$
\end{frame}

\begin{frame}[fragile]{Detailed Balance}
  $$
  \pi\big({\boldsymbol{\theta}}\big)q\big(\boldsymbol{\theta'}, \boldsymbol{\theta}\big) =
  \pi\big({\boldsymbol{\theta'}}\big)q\big(\boldsymbol{\theta}, \boldsymbol{\theta'}\big)
  $$

  Hand waving slightly, we can see that this scheme satisfies the
  premises of the ergodic theorem and so we can conclude that there is
  a unique stationary distribution and $\pi$ must be that
  distribution.
\end{frame}

\subsection{Systematic Scan}

\begin{frame}[fragile]{Systematic Scan}
  Most references on Gibbs sampling do not describe the random scan
  but instead something called a systematic scan.

  Again for simplicity let us consider a model with two parameters. In
  this sampler, we update the parameters in two steps.

  $$
  \begin{aligned}
    \text{Sample } \theta_1^{(i+1)} & \sim \pi(\theta_1 \,\big\vert\, \theta_2^{(i)}) \\
    \text{Sample } \theta_2^{(i+1)} & \sim \pi(\theta_2 \,\big\vert\, \theta_1^{(i+1)})
  \end{aligned}
  $$
\end{frame}

\begin{frame}[fragile]{Systematic Scan}
  \begin{itemize}
    \item
      \textbf{Not} time-homegeneous!!!
      \pause
    \item
    At each step the transition matrix flips between the two
    transition matrices given by the individual steps.
    \pause
  \item
    Cannot apply the ergodic theorem as it only applies to
    time-homogeneous processes.
  \end{itemize}
\end{frame}

\subsection{Bivariate Normal}

\begin{frame}[fragile]{An Example: The Bivariate Normal}
  $$
  \begin{bmatrix}
    \theta_1 \\
    \theta_2
  \end{bmatrix}
  \bigg\vert y \sim
  N \begin{bmatrix}
    \begin{bmatrix}
      y_1 \\
      y_2
    \end{bmatrix} &
    \begin{bmatrix}
      1 & \rho \\
      \rho & 1
    \end{bmatrix}
  \end{bmatrix}
  $$

  The conditional distributions are easily calculated to be

  $$
  \begin{aligned}
    \theta_1 \,\vert\, \theta_2, y &\sim {\cal{N}}(y_1 + \rho(\theta_2 - y_2), 1 - \rho^2) \\
    \theta_2 \,\vert\, \theta_1, y &\sim {\cal{N}}(y_2 + \rho(\theta_1 - y_1), 1 - \rho^2)
  \end{aligned}
  $$
\end{frame}

\begin{frame}[fragile]{Haskell}
\begin{code}
gibbsSampler ::
  Double ->
  RVarT (W.Writer [(Double,Double)]) Double

gibbsSampler oldTheta2 = do
  newTheta1 <- rvarT (Normal (y1 + rho * (oldTheta2 - y2)) var)

  lift $ W.tell [(newTheta1, oldTheta2)]

  newTheta2 <- rvarT (Normal (y2 + rho * (newTheta1 - y1)) var)

  lift $ W.tell [(newTheta1, newTheta2)]

  return $ newTheta2
\end{code}
\end{frame}

\begin{frame}[fragile]{Haskell}
\begin{code}
runMCMC :: Int -> [(Double, Double)]
runMCMC n =
  take n $
  drop burnIn $
  snd $
  runWriter $
  evalStateT (sample (iterateM_ gibbsSampler 2.5))
             (pureMT 2)
\end{code}
\end{frame}

\framedgraphic{Typical Paths}{diagrams/Paths.png}

\framedgraphic{Typical Paths}{diagrams/Samples.png}

\section{Application in Statistics}

\subsection{Known Marginal}

\begin{frame}[fragile]{Prior and Likelihood}
$$
\begin{aligned}
\pi(\mu, \tau) \propto \frac{1}{\tau} & & -\infty < \mu < \infty\, \textrm{and}\, 0 < \tau < \infty
\end{aligned}
$$

$$
p(\boldsymbol{x}\,\vert\,\mu, \sigma) = \prod_{i=1}^n \bigg(\frac{1}{\sigma\sqrt{2\pi}}\bigg)\exp{\bigg( -\frac{(x_i - \mu)^2}{2\sigma^2}\bigg)}
$$
\end{frame}

\begin{frame}[fragile]{Posterior}
Re-writing in terms of precision

$$
p(\boldsymbol{x}\,\vert\,\mu, \tau) \propto \prod_{i=1}^n \sqrt{\tau}\exp{\bigg( -\frac{\tau}{2}{(x_i - \mu)^2}\bigg)}
$$

Thus the posterior is

$$
p(\mu, \tau \,\vert\, \boldsymbol{x}) \propto \tau^{n/2 - 1}\exp{\bigg( -\frac{\tau}{2}\sum_{i=1}^n{(x_i - \mu)^2}\bigg)}
$$
\end{frame}

\begin{frame}[fragile]{Conditionals}
We can re-write the sum in terms of

\begin{itemize}
\uncover<1->{
\item
Sample mean $\bar{x} =
\frac{1}{n}\sum_{i=1}^n x_i$
\item
Sample variance $s^2 =
\frac{1}{n-1}\sum_{i=1}^n (x_i - \bar{x})^2$
\end{itemize}
}
\uncover<2->{
$$
\begin{aligned}
\sum_{i=1}^n (x_i - \mu)^2 &= (n - 1)s^2 + n(\bar{x} - \mu)^2
\end{aligned}
$$
}
\end{frame}

\begin{frame}[fragile]{Conditionals}
Thus the conditional posterior for $\mu$ is

$$
\begin{aligned}
p(\mu \,\vert\, \tau, \boldsymbol{x}) &\propto \exp{\bigg( -\frac{\tau}{2}\bigg(\nu s^2 + \sum_{i=1}^n{(\mu - \bar{x})^2}\bigg)\bigg)} \\
&\propto \exp{\bigg( -\frac{n\tau}{2}{(\mu - \bar{x})^2}\bigg)} \\
\end{aligned}
$$

which we recognise as a normal distribution with mean of $\bar{x}$ and
a variance of $(n\tau)^{-1}$.
\end{frame}

\begin{frame}[fragile]{Conditionals}
The conditional posterior for $\tau$ is

$$
\begin{aligned}
p(\tau \,\vert\, , \mu, \boldsymbol{x}) &\propto \tau^{n/2 -1}\exp\bigg(-\tau\frac{1}{2}\sum_{i=1}^n{(x_i - \mu)^2}\bigg)
\end{aligned}
$$

which we recognise as a gamma distribution with a shape of $n/2$ and a scale of $\frac{1}{2}\sum_{i=1}^n{(x_i - \mu)^2}$
\end{frame}

\begin{frame}[fragile]{Marginal}
In this particular case, we can calculate the marginal posterior of
$\mu$ analytically.

$$
\begin{aligned}
p(\mu \,\vert\, \boldsymbol{x})
&\propto \bigg( 1 + \frac{n(\mu - \bar{x})^2}{(n - 1)s^2} \bigg)^{-n/2} \\
\end{aligned}
$$

This is the
non-standardized Student's t-distribution $t_{n-1}(\bar{x}, s^2/n)$.
\end{frame}

\begin{frame}[fragile]{Haskell}
\begin{code}
gibbsSampler :: MonadRandom m =>
                Double -> m (Maybe ((Double, Double), Double))
gibbsSampler oldTau = do
  newMu <- sample (Normal xBar (recip (sqrt (n * oldTau))))
  let shape = 0.5 * n
      scale = 0.5 * (x2Sum + n * newMu^2 - 2 * n * newMu * xBar)
  newTau <- sample (Gamma shape (recip scale))
  return $ Just ((newMu, newTau), newTau)
\end{code}

From which we can create an infinite stream of samples.
\begin{code}
gibbsSamples :: [(Double, Double)]
gibbsSamples = evalState (unfoldrM gibbsSampler initTau)
                         (pureMT 1)
\end{code}
\end{frame}

\framedgraphic{Posterior for $\mu$}{diagrams/Posterior.png}

\begin{frame}[fragile]{JAGS}
JAGS is a mature, DSL for building Bayesian statistical models
using Gibbs(?) sampling.

Here is our model as expressed in JAGS. Somewhat terse.
\begin{scriptsize}
\begin{lstlisting}
model {
	for (i in 1:N) {
		x[i] ~ dnorm(mu, tau)
	}
	mu ~ dnorm(0, 1.0E-6)
	tau <- pow(sigma, -2)
	sigma ~ dunif(0, 1000)
}
\end{lstlisting}
\end{scriptsize}
\end{frame}

\framedgraphic{Output for JAGS}{diagrams/jags.png}

\begin{frame}[fragile]{STAN} STAN: DSL similar to JAGS but newer, uses
HMC, allows variable re-assignment, cannot really be described as
declarative.
\begin{scriptsize}
\begin{lstlisting}
data {
  int<lower=0> N;
  real x[N];
}
parameters {
  real mu;
  real<lower=0,upper=1000> sigma;
}
model{
  x     ~ normal(mu, sigma);
  mu    ~ normal(0, 1000);
}
\end{lstlisting}
\end{scriptsize}
\end{frame}

\framedgraphic{Output for JAGS}{diagrams/stan.png}

\begin{frame}[fragile]{Resources}
\begin{itemize}
\item \textcolor{blue}{http://idontgetoutmuch.wordpress.com}
\item
  \textcolor{blue}{http://en.wikipedia.org/wiki/Automatic\_differentiation}
\item \textcolor{blue}{http://www.autodiff.org}
\item \textcolor{blue}{http://hackage.haskell.org/package/ad}
\item \textcolor{blue}{http://mc-stan.org}
\item \textcolor{blue}{The Monad.Reader Issue 21}
\end{itemize}
\end{frame}

\end{document}
