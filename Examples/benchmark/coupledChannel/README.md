# Coupled benchmark analysis

a) For perfoming this benchmark analysis, download the relevant toy data
samples:

1. **pbar p -> eta eta pi0** <br> data events at: <br>
   http://www.ep1.rub.de/~bertram/benchmarkCoupledChannel/genPbarpToEtaEtaPi0Data.dat
   <br> Monte Carlo events at: <br>
   http://www.ep1.rub.de/~bertram/benchmarkCoupledChannel/genPbarpToEtaEtaPi0Mc.dat

2. **pbar p -> K+ K- pi0** <br> data events at: <br>
   http://www.ep1.rub.de/~bertram/benchmarkCoupledChannel/evtGenPbarpToKKPiData.dat
   <br> Monte Carlo events at: <br>
   http://www.ep1.rub.de/~bertram/benchmarkCoupledChannel/evtGenPbarpToKKPiMc.dat

3. **e+ e- -> J/psi -> chi_c1 gamma -> (K+ K- pi0) gamma** <br> data events at:
   <br>
   http://www.ep1.rub.de/~bertram/benchmarkCoupledChannel/evtGenChic1ToKKpiData.dat
   <br> Monte Carlo events at: <br>
   http://www.ep1.rub.de/~bertram/benchmarkCoupledChannel/evtGenChic1ToKKpiMc.dat

4. **D0 -> K- pi+ pi0** <br> data events at: <br>
   http://www.ep1.rub.de/~bertram/benchmarkCoupledChannel/evtGenD0ToKmPipPi0Data.dat
   <br> Monte Carlo events at: <br>
   http://www.ep1.rub.de/~bertram/benchmarkCoupledChannel/evtGenD0ToKmPipPi0Mc.dat

On linux machines, you can run the script
[`downloadData_script`](./downloadData_script) for downloading all required
event files.

b) By starting the script
[`benchmarkCoupledChannel_runscript`](./benchmarkCoupledChannel_runscript), the
PWA will automatically be performed in the server-client mode with 4 clients
(with one client per channel).

c) Compare the result with the reference log file `logOutput.server`.

d) Start `"qaMode"` for the individual channels with
[`qa_runscript`](./qa_runscript).

e) Compare data with fit result by starting ROOT:

```
root [0] .L \$TOP_DIR/Scripts/DataFitCompare.C+g
root [1] DataFitCompareAll("pawianHists_pbarpEtaEtaPi.root")
root [2] DataFitCompareAll("pawianHists_pbarpKKPi.root")
root [3] DataFitCompareAll("pawianHists_JpsiToChic1GammaToKKPiGamma.root")
root [4] DataFitCompareAll("pawianHists_D0ToKPiPi.root")
```
