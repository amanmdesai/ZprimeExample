import ROOT
import scipy.stats as sp
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.special import factorial

f = ROOT.TFile.Open("DataSample_search.root")

h_mass_data = f.Get('h_mass_data')
h_mass_bgr = f.Get('h_mass_bgr')
h_mass_sig_100 = f.Get('h_mass_sig_100')
h_mass_sig_250 = f.Get('h_mass_sig_250')
h_mass_sig_400 = f.Get('h_mass_sig_400')

h_mass_bgr.SetFillColor(ROOT.kYellow)
h_mass_sig_100.SetFillColor(ROOT.kBlue)
h_mass_sig_250.SetFillColor(ROOT.kRed)
h_mass_sig_400.SetFillColor(ROOT.kCyan)

#h_mass_bgr.SetAxisRange(0, 500, "X")
#h_mass_sig_100.SetAxisRange(0, 500, "X")
#h_mass_sig_250.SetAxisRange(0, 500, "X")
#h_mass_sig_400.SetAxisRange(0, 500, "X")
#_mass_data.SetAxisRange(0, 500, "X")
c1 = ROOT.TCanvas()
stack = ROOT.THStack()

"""
h_mass_bgr.GetXaxis().SetRangeUser(0, 500)
h_mass_sig_100.GetXaxis().SetRangeUser(0, 500)
h_mass_sig_250.GetXaxis().SetRangeUser(0, 500)
h_mass_sig_400.GetXaxis().SetRangeUser(0, 500)
h_mass_data.GetXaxis().SetRangeUser(0, 500)
"""

h_mass_data.SetMarkerStyle(20)
h_mass_data.SetMarkerColor(ROOT.kBlack)


stack.Add(h_mass_bgr)
stack.Add(h_mass_sig_100)
stack.Add(h_mass_sig_250)
stack.Add(h_mass_sig_400)

#h_mass_sig_100.Draw("hist")
#h_mass_sig_250.Draw("histsame")
#h_mass_sig_400.Draw("histsame")
#stack.Draw("histsame")
h_mass_data.Sumw2()
stack.Draw("hist")
h_mass_data.Draw("e1same")
#c1.Update()
#ROOT.gPad.Modified()
#ROOT.gPad.Update()
c1.SaveAs("zprime.pdf")
mass_bin = 50.
exp_bkg = h_mass_bgr.Integral(h_mass_bgr.FindBin(250-.5*mass_bin),h_mass_bgr.FindBin(250+.5*mass_bin))
exp_signal = h_mass_sig_250.Integral(h_mass_sig_250.FindBin(250-.5*mass_bin),h_mass_sig_250.FindBin(250+.5*mass_bin))
obs_data = h_mass_data.Integral(h_mass_data.FindBin(250-.5*mass_bin),h_mass_data.FindBin(250+.5*mass_bin))

### significance

def significance(bkg, sig, data):
    Z = 2*data*ROOT.TMath.Log(1+(sig/bkg)) - 2*sig
    Z = ROOT.TMath.Sqrt(Z)
    return Z


#generate figure 11.2 b
plt.hist(sp.poisson.rvs(mu=14.71, size=1000000),histtype="step",bins=25)
plt.hist(sp.poisson.rvs(mu=22.55, size=1000000),histtype="step",bins=25)
#plt.show()


### side band fit
value_1 = 150.
value_2 = 350.


h_mass_data_rebin = f.Get('h_mass_data').Clone("h_mass_data_rebin")
h_mass_bgr_rebin = f.Get('h_mass_bgr').Clone("h_mass_bgr_rebin")
h_mass_sig_250_rebin = f.Get('h_mass_sig_250').Clone("h_mass_sig_250_rebin")

rebin_factor = 10
h_mass_bgr_rebin.Rebin(rebin_factor)
h_mass_sig_250_rebin.Rebin(rebin_factor)
h_mass_data_rebin.Rebin(rebin_factor)

Num_bins = 2*mass_bin/h_mass_bgr_rebin.GetBinWidth(1)

bin_bgr = np.zeros(int(Num_bins)*2)
bin_sig = np.zeros(int(Num_bins)*2)
bin_data = np.zeros(int(Num_bins)*2)

for i in range(int(Num_bins)):
    bin_bgr[i] = h_mass_bgr_rebin.GetBinContent(h_mass_bgr_rebin.FindBin(value_1-mass_bin)+i)
    bin_sig[i] = h_mass_sig_250_rebin.GetBinContent(h_mass_sig_250_rebin.FindBin(value_1-mass_bin)+i)
    bin_data[i] = h_mass_data_rebin.GetBinContent(h_mass_data_rebin.FindBin(value_1-mass_bin)+i)


for i in range(int(Num_bins)-1,int(Num_bins)*2):
    bin_bgr[i] = h_mass_bgr_rebin.GetBinContent(h_mass_bgr_rebin.FindBin(value_2-mass_bin)+i)
    bin_sig[i] = h_mass_sig_250_rebin.GetBinContent(h_mass_sig_250_rebin.FindBin(value_2-mass_bin)+i)
    bin_data[i] = h_mass_data_rebin.GetBinContent(h_mass_data_rebin.FindBin(value_2-mass_bin)+i)

def _Log(alpha, bi, ni):
    result=0
    for i in range(len(bi)):
        result += - ni[i]*np.log(alpha*bi[i]) + bi[i]*alpha + np.log(factorial(ni[i]))
    return result*2

num_alpha = 200
alpha = np.linspace(0.7,1.4,num_alpha)
side_band_fit = np.zeros(num_alpha)

for i in range(num_alpha):
    side_band_fit[i] = _Log(alpha[i], bin_bgr, bin_data)

idx = np.argwhere(np.diff(np.sign(side_band_fit - np.ones(len(side_band_fit))*(min(side_band_fit)+1)))).flatten()

plt.plot(alpha, side_band_fit)
#plt.plot(alpha, np.ones(len(side_band_fit))*(min(side_band_fit)+1))
plt.plot(alpha[idx],side_band_fit[idx], marker='o')
plt.xlim(0.8,1.3)
plt.ylim(side_band_fit.min()-3,side_band_fit.min()+3)
#plt.show()


#print(np.argmin(side_band_fit+1))
#print(alpha[np.argmin(side_band_fit+1)])
#print(side_band_fit -  np.ones(len(side_band_fit))*(min(side_band_fit)+1))

##region 1
side_1_bgr = h_mass_bgr.Integral(h_mass_bgr.FindBin(value_1-mass_bin),h_mass_bgr.FindBin(value_1+mass_bin))
side_1_data = h_mass_data.Integral(h_mass_data.FindBin(value_1-mass_bin),h_mass_data.FindBin(value_1+mass_bin))
##region 2
side_2_bgr = h_mass_bgr.Integral(h_mass_bgr.FindBin(value_2-mass_bin),h_mass_bgr.FindBin(value_2+mass_bin))
side_2_data = h_mass_data.Integral(h_mass_data.FindBin(value_2-mass_bin),h_mass_data.FindBin(value_2+mass_bin))
calc_alpha =(side_1_data+side_2_data)/(side_1_bgr+side_2_bgr)
new_exp_bkg = exp_bkg*alpha[np.argmin(side_band_fit)]

### side band ends here


### t statistics evaluation begins here

value_3 = 350
mass_bin2 = 250
#Num_bins = 2*mass_bin/h_mass_bgr_rebin.GetBinWidth(1)
Num_bins = (h_mass_data_rebin.FindBin(600) - h_mass_data_rebin.FindBin(100)) #2*mass_bin2/h_mass_bgr_rebin.GetBinWidth(1)

bin_bgr_full = np.zeros(int(Num_bins))
bin_sig_full = np.zeros(int(Num_bins))
bin_data_full = np.zeros(int(Num_bins))
#print(h_mass_data_rebin.FindBin(100),h_mass_data_rebin.FindBin(600), Num_bins)
for i in range(int(Num_bins)):
    bin_bgr_full[i] = h_mass_bgr_rebin.GetBinContent(h_mass_bgr_rebin.FindBin(value_3-mass_bin2)+i)
    bin_sig_full[i] = h_mass_sig_250_rebin.GetBinContent(h_mass_sig_250_rebin.FindBin(value_3-mass_bin2)+i)
    bin_data_full[i] = h_mass_data_rebin.GetBinContent(h_mass_data_rebin.FindBin(value_3-mass_bin2)+i)

#plt.hist(bin_data_full,histtype='step')#,linestyle='none',marker='o')
#plt.hist(bin_bgr_full,histtype='step')
#plt.hist(bin_sig_full+bin_bgr_full,histtype='step')
##plt.show()

def _Log2(bi,ni):
    result=0
    for i in range(len(ni)):
        #        result += - ni[i]*np.log(alpha*bi[i] + mu*si[i]) + (bi[i]*alpha + mu*si[i]) + np.log(factorial(ni[i]))
        result += - 2.*ni[i]*np.log(bi[i]) + 2*bi[i] + 2*np.log(factorial(ni[i]))
    return result

num_alpha = 100
alpha2 = np.linspace(0.01,5.,num_alpha)
mu = np.linspace(0.0,5.,num_alpha)
_fit_0 = np.zeros(num_alpha)
_fit_1 = np.zeros(num_alpha)
_full_fit = np.zeros((num_alpha,num_alpha))

#X, Y = np.meshgrid(alpha2, mu)
##attempt to make 2d contour for Log2 with variables mu and alpha
#for i in range(num_alpha):
#    for j in range(num_alpha):
#        _full_fit[i][j] = _Log2(alpha[i],mu[j], bin_bgr_full,bin_sig_full, bin_data_full)
# plt.contour(X, Y,_full_fit)

for i in range(num_alpha):
    _fit_0[i] = _Log2(alpha2[i]*bin_bgr_full, bin_data_full)
    _fit_1[i] = _Log2(alpha2[i]*bin_bgr_full + bin_sig_full, bin_data_full)

t = -_fit_1 + _fit_0
alpha_val = alpha[np.argmin(side_band_fit)]

_fit_0 = _Log2(bin_bgr_full, bin_data_full)
_fit_1 = _Log2(bin_bgr_full + bin_sig_full, bin_data_full)

t_obs = -_fit_1 + _fit_0

#plt.plot(alpha2,_fit_0)
#plt.plot(alpha2,_fit_0)
#plt.plot(alpha2,_fit_1)
#plt.plot(alpha2,t)
#plt.xlim(0.1,3.2)
##plt.show()
#plt.plot(alpha2,_fit_0)
#plt.plot(alpha2,_fit_1)
#plt.ylim(side_band_fit.min()-5,side_band_fit.min()+5)
#plt.plot(t,_fit_1)
#plt.hist(t,bins=100)
##plt.show()
#plt.xlim(0,2)
#plt.ylim(0,2)

### t statistics evaluation ends here

### confidence limits

#print(t)

# First Draw events from the main histogram and generate poisson random distribution, evaluate the t statistics for this data

array_data_b = np.zeros(len(bin_bgr_full))
array_data_sb = np.zeros(len(bin_bgr_full))



def calc_toy():
    for i in range(len(bin_bgr_full)):
        #mass_data.SetBinContent(i,np.random.poisson(bgr_full.GetBinContent(i)))
        #mass_data2.SetBinContent(i,np.random.poisson(bgr_full.GetBinContent(i)+sig_full.GetBinContent(i)))

        array_data_b[i] = np.random.poisson(bin_bgr_full[i])#.GetBinContent(i))
        array_data_sb[i] = np.random.poisson(bin_bgr_full[i]+bin_sig_full[i])#.GetBinContent(i))
        #array_data_sb[i-1] = np.random.poisson(bgr_full.GetBinContent(i)+sig_full.GetBinContent(i))
    return array_data_b, array_data_sb

def t_stat():
    array_data_b, array_data_sb = calc_toy()
    _fit_0_b=  _Log2(bin_bgr_full, array_data_b)
    _fit_1_b = _Log2(bin_bgr_full+ bin_sig_full, array_data_b)
    _fit_0_s = _Log2(bin_bgr_full, array_data_sb)
    _fit_1_s = _Log2(bin_bgr_full+ bin_sig_full, array_data_sb)
    t_b = +_fit_0_b - _fit_1_b
    t_s = +_fit_0_s - _fit_1_s
    return t_b, t_s

size = 5000
tval_s = np.zeros(size)
tval_b = np.zeros(size)
for i in range(size):
    tval_s[i], tval_b[i] = t_stat()


plt.show()

#print(.2*np.sum(abs(tval_b[tval_b > t_obs]))/size)
#print(.2*np.sum(abs(tval_s[tval_s < t_obs]))/size)


plt.xlim(-15,25)
ns ,sbins, _ = plt.hist(tval_s,alpha = .6,label='s+b',bins=30,density =True)
nb ,bbins, _ = plt.hist(tval_b,alpha = .6,label='b-only',bins=30,density =True)
plt.plot(t_obs,.0,marker='x',linestyle='None',color='Black',markersize=20,label='Observed')
plt.plot(np.median(tval_b),.0,marker='o',linestyle='None',color='Black',markersize=20,label='Observed')
plt.plot(np.median(tval_s),.0,marker='o',linestyle='None',color='Black',markersize=20,label='Observed')
plt.legend()
plt.show()

print(sbins[0],bbins[0])

bin_int_s = int((np.median(tval_s) + sbins[0]) /(30/40))
bin_int_b = int((np.median(tval_b) + bbins[0]) /(30/40))
bin_int_data = int((np.median(t_obs) +15) /(30/40))


bin_width = sbins[1] - sbins[0]
# sum over number in each bin and mult by bin width, which can be factored out

integral_CLb_1 = bin_width * sum(nb[:bin_int_s])#/sum(ns)
integral_CLsb_1 = bin_width * sum(ns[bin_int_s:])#/sum(ns)


integral_CLb_2 = bin_width * sum(nb[:bin_int_b])#/sum(nb)
integral_CLsb_2 = bin_width * sum(ns[bin_int_b:])#/sum(nb)

integral_CLb_data = bin_width * sum(nb[:bin_int_data])#/sum(nb)
integral_CLsb_data = bin_width * sum(ns[bin_int_data:])#/sum(nb)


#new_bgr_sig[i] = np.random.poisson(bgr_full[i-1]+sig_full[i-1])
#plt.hist(sp.poisson.rvs(mu=exp_bkg, size=1000000),histtype="step",bins=25)
#plt.hist(sp.poisson.rvs(exp_bkg+exp_signal, size=1000000),histtype="step",bins=25)
##plt.show()

#print(1 - sp.poisson.cdf(obs_data-1,exp_bkg+exp_signal))



####          RESULTS

print("|====================================================================|")
print("|------------------  Analysis Report Begins -------------------------|")
print("|====================================================================|")
print("|------------------  Sample Statistics ------------------------------|")
print("| + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +|")
print("| Expected number of background events     | ",exp_bkg )
print("|--------------------------------------------------------------------|")
print("| Expected number of signal events         | ",exp_signal)
print("|--------------------------------------------------------------------|")
print("| Sum of Expected Signal+background events | ",exp_bkg+exp_signal )
print("|--------------------------------------------------------------------|")
print("| Observed Data                            | ",obs_data)
print("|====================================================================|")
print("|-----------------------   P VALUE  ---------------------------------|")
print("| + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +|")
print("| Expected p value                         | ",1 - sp.poisson.cdf(exp_bkg+exp_signal-1,exp_bkg))
print("|--------------------------------------------------------------------|")
print("| Observed p value                         | ",1 - sp.poisson.cdf(obs_data-1,exp_bkg))
print("|====================================================================|")
print("|-----------------------Significance --------------------------------|")
print("| + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +|")
print("| Expected significance                    | ",significance(exp_bkg, exp_signal, exp_bkg+exp_signal ))
print("|--------------------------------------------------------------------|")
print("| Observed significance                    | ",significance(exp_bkg, exp_signal, obs_data))
print("|====================================================================|")
print("|-----------------------Side Band Fit -------------------------------|")
print("| + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +|")
print("| Calculation                              | ", alpha[np.argmin(side_band_fit)]," - ",  alpha[idx[0]]-alpha[np.argmin(side_band_fit)], " + ", alpha[idx[1]]-alpha[np.argmin(side_band_fit)]) # (alpha[np.argmin(side_band_fit)]/math.sqrt(len(bin_bgr))))
print("|--------------------------------------------------------------------|")
print("| Analytical calculation                   | ", calc_alpha,    " +- ", calc_alpha/math.sqrt(len(bin_bgr)))
print("|--------------------------------------------------------------------|")
print("| Expected Background in Signal region     | ",  new_exp_bkg)
print("|--------------------------------------------------------------------|")
print("| New Observed significance                | ",significance(new_exp_bkg, exp_signal, obs_data))
print("|--------------------------------------------------------------------|")
print("| New Expected significance                | ",significance(new_exp_bkg, exp_signal, new_exp_bkg+exp_signal ))
print("|====================================================================|")
print("|------------------  Confidence Limits ------------------------------|")
print("| Observed t value                         | ",t_obs)
print("|--------------------------------------------------------------------|")
print("| Median t value for S+B                   | ",np.median(tval_s))
print("|--------------------------------------------------------------------|")
print("| Median t value for B                     | ",np.median(tval_b))
print("|--------------------------------------------------------------------|")
print("| 1 - CL b (Median s+b)                    | ",integral_CLb_1)
print("|--------------------------------------------------------------------|")
print("| CL s+b (Median s+b)                      | ",integral_CLsb_1)
print("|--------------------------------------------------------------------|")
print("| 1 - CL b (Median b)                      | ",integral_CLb_2)
print("|--------------------------------------------------------------------|")
print("| CL s+b (Median b)                        | ",integral_CLsb_2)
print("|--------------------------------------------------------------------|")
print("| 1 - CL b (Data)                          | ",integral_CLb_data)
print("|--------------------------------------------------------------------|")
print("| CL s+b (Data)                            | ",integral_CLsb_data)
print("|--------------------------------------------------------------------|")
print("|====================================================================|")
print("|------------------  Analysis Report Ends ---------------------------|")
print("|====================================================================|")


## Numbers from the book
"""
print("expected significance ",1 - sp.poisson.cdf(22.55, 14.71, loc=0)) # expected because we used sum of signal + background
print("observed significance ",1 - sp.poisson.cdf(25-1, 14.71, loc=0)) # observed because we used 25 data events

def test_poisson(n, v):
    result = pow(v,n)*math.exp(-v)
    result = result/math.factorial(int(n))
    return result

def test_Log(alpha, bi, ni):
    result =0.
    for i in range(len(bi)):
        result += poisson(ni[i], alpha*bi[i])
    return result

"""
