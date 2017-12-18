import numpy as np
import matplotlib.pyplot as plt

n = 100000
k = 100
guess_0 = np.zeros(n)
guess = np.zeros((k,n))
for i in np.arange(k):
    guess[i,::i+1] = 1.
answers = np.zeros(n)
answers[::8] = 1.
np.random.shuffle(answers)

score_0 = np.sum(np.equal(guess_0,answers))
score_0 = score_0/float(n)
score = np.zeros(k)
for i in np.arange(k):
    score[i] = np.sum(np.equal(guess[i, :],answers))
    score[i] = score[i]/float(n)
print 'score_0 = ' + str(score_0)
for i in np.arange(k):
    print 'score ' + str(i) + ' = ' + str(score[i])
plt.plot(score)
plt.axhline(score_0, ls='--')
plt.show()
