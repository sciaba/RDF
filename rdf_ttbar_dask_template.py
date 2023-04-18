NUM_WORKERS = _WORKERS_
DATASETS = "_DATASETS_"

import os
import re
import json
import time
from urllib.request import urlretrieve

from ROOT import TCanvas, THStack
import ROOT

from distributed import Client, LocalCluster

RDataFrame = ROOT.RDF.Experimental.Distributed.Dask.RDataFrame
RunGraphs = ROOT.RDF.Experimental.Distributed.RunGraphs
VariationsFor = ROOT.RDF.Experimental.Distributed.VariationsFor
initialize = ROOT.RDF.Experimental.Distributed.initialize

def create_connection() -> Client:
    client = Client(
        LocalCluster(n_workers=NUM_WORKERS, threads_per_worker=1,
                     processes=True)
    )
    return client

def myinit():
    ROOT.gSystem.CompileMacro("helper.cpp", "kO")
    ROOT.gInterpreter.Declare(f"""
    #ifndef MYPTR
    #define MYPTR
    auto pt_res_up_obj = pt_res_up(1);
    #endif
    """)

print(f'The num of workers = {NUM_WORKERS}')

N_FILES_MAX_PER_SAMPLE = _NFILES_
DOWNLOAD = False
LOCAL = False
AF_NAME = "_AFNAME_"
FILE = f'rdf-{N_FILES_MAX_PER_SAMPLE}.root'

class TtbarAnalysis(dict):

    def __init__(self, n_files_max_per_sample=1, num_bins=25, bin_low=50,
                 bin_high=550, download_input_data=False,
                 use_local_data=False, af_name="unl", datasets="ntuples.json"):
        
        self.variations = {}
        self.download_input_data = download_input_data
        self.use_local_data = use_local_data
        self._nevts_total = {}
        self.n_files_max_per_sample = n_files_max_per_sample
        self.input_data = self._construct_fileset(af_name, datasets,
                                                  use_xcache=_XCACHE_)
        self.num_bins = num_bins
        self.bin_low = bin_low
        self.bin_high = bin_high
        self.xsec_info = {
            "ttbar": 396.87 + 332.97,
            "single_top_s_chan": 2.0268 + 1.2676,
            "single_top_t_chan": (36.993 + 22.175)/0.252,
            "single_top_tW": 37.936 + 37.906,
            "wjets": 61457 * 0.252,
            "data": None
        }

    def _construct_fileset(self, af_name, datasets, use_xcache=False):
        n_files_max_per_sample = self.n_files_max_per_sample
        with open (datasets) as f:
            file_info = json.load(f)
        fileset = {}
        for process in file_info.keys():
            if process == "data":
                continue
            fileset[process] = {}
            self[process] = {}
            self._nevts_total[process] = {}
            for variation in file_info[process].keys():
                file_list = file_info[process][variation]["files"]
                if n_files_max_per_sample != -1:
                    file_list = file_list[:n_files_max_per_sample]
                file_paths = [f["path"] for f in file_list]
                if af_name == "unl-xrootd":
                    orig_prefix = "https://xrootd-local.unl.edu:1094"
                    new_prefix = "root://xrootd-local.unl.edu:1094"
                    file_paths = [f.replace(orig_prefix,
                                            new_prefix) for f in file_paths]
                orig_prefix = "https://xrootd-local.unl.edu:1094//store/user/AGC"
                if (af_name == "cern-xrootd"):
                    if (re.search("merged", datasets)):
                        new_prefix = "root://eoscms.cern.ch//eos/cms/store/test/agc"
                        file_paths = [f.replace(orig_prefix,
                                                new_prefix) for f in file_paths]
                    else:
                        new_prefix = "root://eoscms.cern.ch//eos/cms/opstest/asciaba/agc"
                        file_paths = [f.replace(orig_prefix,
                                                new_prefix) for f in file_paths]
                if (af_name == "cernbox-xrootd"):
                    new_prefix = "root://eosuser.cern.ch//eos/user/a/asciaba/datasets/agc"
                    file_paths = [f.replace(orig_prefix,
                                            new_prefix) for f in file_paths]
                if (af_name == "cern-local"):
                    new_prefix = "/data/datasets/agc"
                    file_paths = [f.replace(orig_prefix,
                                            new_prefix) for f in file_paths]
                if use_xcache:
                    file_paths = [f.replace("root:", "root://xcache01.cern.ch//xroot:") for f in file_paths]
                nevts_total = sum([f["nevts"] for f in file_list])
                self._nevts_total[process].update({variation:nevts_total})
                fileset[process].update({variation: file_paths})
                if (self.download_input_data):
                    dir_name = f"input/{process}_{variation}"
                    os.makedirs(dir_name, exist_ok=True)
                    for i in range(len(file_paths)):
                        path = file_paths[i]
                        file = f"{dir_name}/{i}.root"
                        if not os.path.exists(file):
                            urlretrieve(path, file)   
                            print(f"{file} has been created")
                        else:
                            print(f"{file} already exists")
                self[process][variation] = {}
                
        return fileset

    def fill(self, process, variation, connection):
        
        input_data = self.input_data[process][variation]               
        d = RDataFrame('events', input_data, daskclient=connection,
                       npartitions=NUM_WORKERS)
        d._headnode.backend.distribute_unique_paths(["helper.cpp"])
        
        # normalization for MC
        x_sec = self.xsec_info[process]
        nevts_total = self._nevts_total[process][variation]
        lumi = 3378
        xsec_weight = x_sec * lumi / nevts_total
        d = d.Define('weights', str(xsec_weight))
        
        if variation == 'nominal':            
            d = d.Vary('jet_pt', "ROOT::RVec<ROOT::RVecF>{jet_pt*pt_scale_up(), jet_pt*pt_res_up_obj(jet_pt, rdfslot_)}", ["pt_scale_up", "pt_res_up"])
            if process == 'wjets':                
                d = d.Vary('weights', 
                           "weights*flat_variation()",
                           [f"scale_var_{direction}" for direction in ["up", "down"]]
                          )
                
        d = d.Define('electron_pt_mask', 'electron_pt>25').Define('muon_pt_mask', 'muon_pt>25').Define('jet_pt_mask', 'jet_pt>25')\
             .Filter('Sum(electron_pt_mask) + Sum(muon_pt_mask) == 1')\
             .Filter('Sum(jet_pt_mask) >= 4')\
             .Filter('Sum(jet_btag[jet_pt_mask]>=0.5)>=1')
             
        
        d = d.Vary('weights', 
                   'ROOT::RVecD{weights*btag_weight_variation(jet_pt[jet_pt_mask])}',
                   [f"{weight_name}_{direction}" for weight_name in [f"btag_var_{i}" for i in range(4)] for direction in ["up", "down"]]
                  ) if variation == 'nominal' else d
        
        measured = {"4j1b": "HT", "4j2b": 'trijet_mass'} # columns names of observables for two regions
        for region in ["4j1b","4j2b"]:
            observable = measured[region]
            
            if region == "4j1b":
                
                fork = d.Filter('Sum(jet_btag[jet_pt_mask]>=0.5)==1').Define(observable, 'Sum(jet_pt[jet_pt_mask])')      

            elif region == "4j2b":
                
                fork = d.Filter('Sum(jet_btag[jet_pt_mask]>=0.5)>1').Define("jet_p4", 
                    "ROOT::VecOps::Construct<ROOT::Math::PxPyPzMVector>(jet_px[jet_pt_mask], jet_py[jet_pt_mask], jet_pz[jet_pt_mask], jet_mass[jet_pt_mask])"
                )
                
                fork = fork.Define('trijet', 
                    'ROOT::VecOps::Combinations(jet_pt[jet_pt_mask],3)'
                ).Define('ntrijet', 'trijet[0].size()')

                fork = fork.Define('trijet_p4', 
                                      'ROOT::VecOps::RVec<ROOT::Math::PxPyPzMVector> trijet_p4(ntrijet);'              +\
                                      'for (int i = 0; i < ntrijet; ++i) {'                                            +\
                                          'int j1 = trijet[0][i]; int j2 = trijet[1][i]; int j3 = trijet[2][i];'       +\
                                          'trijet_p4[i] = jet_p4[j1] + jet_p4[j2] + jet_p4[j3];'                       +\
                                      '}'                                                                              +\
                                      'return trijet_p4;'                                                                                                                          
                                     )

                fork = fork.Define('trijet_pt', 
                        'return ROOT::VecOps::Map(trijet_p4, [](ROOT::Math::PxPyPzMVector v) { return v.Pt(); })'
                                            )

                fork = fork.Define('trijet_btag', 
                                                  'ROOT::VecOps::RVec<bool> btag(ntrijet);'                                   +\
                                                  'for (int i = 0; i < ntrijet; ++i) {'                                       +\
                                                   'int j1 = trijet[0][i]; int j2 = trijet[1][i]; int j3 = trijet[2][i];'     +\
                                                   'btag[i]=std::max({jet_btag[j1], jet_btag[j2], jet_btag[j3]})>0.5;'        +\
                                                  '}'                                                                         +\
                                                  'return btag;'
                                            )
                fork=fork.Define(observable,
                                                  'double mass;'+\
                                                  'double Pt = 0;'+\
                                                  'double indx = 0;'+\
                                                  'for (int i = 0; i < ntrijet; ++i) {'               +\
                                                  '    if ((Pt < trijet_pt[i]) && (trijet_btag[i])) {'+\
                                                  '        Pt = trijet_pt[i];'+\
                                                  '        indx=i;'+\
                                                  '    }'                                            +\
                                                  '}'                                                +\
                                                  'mass = trijet_p4[indx].M();'             +\
                                                  'return mass;'
                                                 )
                
            res = fork.Histo1D((f'{process}_{variation}_{region}', process,
                                self.num_bins, self.bin_low, self.bin_high), observable, 'weights')
            self.hist.append(res)
            print(f'histogram {region}_{process}_{variation} has been created')
            
            if variation == 'nominal':
                self.variations[f"{process}__{region}"] = VariationsFor(res)
            else:
                self[process][variation][region] = res

    def Fill(self, connection):
        self.hist = []
        for process in self:
            
            for variation in self.input_data[process]:
                self.fill(process=process, variation=variation,
                          connection=connection)

    def Accumulate(self):
        RunGraphs(self.hist)  
    
    def TransfToDict(self):
        for key in self.variations.keys():
            hist_map = self.variations[key]
            key = str(key).split('__')
            process = key[0]; region = key[1]
            for hist_name in hist_map.GetKeys():
                variation = 'nominal' if hist_name == 'nominal' else str(hist_name).split(':')[1]
                if variation not in self[process]:
                    self[process][variation] = {}
                hist = hist_map[hist_name]
                if not isinstance(hist, ROOT.TH1D):
                    hist = hist.GetValue()
                self[process][variation][region] = hist
        self.ExportJSON()
        
    def GetProcStack(self, region, variation='nominal'):
        return [self[process][variation][region] for process in self]
    
    def GetVarStack(self, region, process="ttbar"):
        ret = [self[process][variation][region] for variation in self[process]]
        ret = [h.GetValue() for h in ret if not isinstance(h, ROOT.TH1D)]
        return ret
        
    def ExportJSON(self):
        data = {}
        for process in self:
            data[process] = {}
            for variation in self[process]:
                data[process][variation] = [region for region in self[process][variation]]
        with open('data.json', 'w') as f:
            json.dump(data, f)
                
def analyse():
    initialize(myinit)
    analysisManager = TtbarAnalysis(download_input_data=DOWNLOAD,
                                    use_local_data=LOCAL,
                                    n_files_max_per_sample = N_FILES_MAX_PER_SAMPLE,
                                    af_name = AF_NAME,
                                    datasets=DATASETS)

    print(f"processes in fileset: {list(analysisManager.keys())}")
    print(f"\nexample of information inside analysisManager:\n{{\n  'urls': [{analysisManager.input_data['ttbar']['nominal'][0]}, ...],")

    with create_connection() as conn:
        t0 = time.time()
        analysisManager.Fill(connection=conn)
        t1 = time.time()
        print(f"\npreprocessing took {round(t1 - t0,2)} seconds")
        analysisManager.Accumulate()
        t2 = time.time()
        print(f"processing took {round(t2 - t1,2)} seconds")
        print(f"execution took {round(t2 - t0,2)} seconds")

    analysisManager.TransfToDict()
    return analysisManager

def make_plots(analysisManager):
    c = TCanvas('c', 'c', 3000*2, 2000*2) 
    hlist = analysisManager.GetProcStack(region='4j1b')
    hs = THStack('j4b1', '>=4 jets, 1 b-tag; H_{T} [GeV]')
    for h in hlist:
        h = ROOT.Slice(h, 120, 550)
        ptr = h.Rebin(2, h.GetTitle())
        hs.Add(ptr)
    hs.Draw('hist pfc plc')
    c.Draw()
    x = hs.GetXaxis()
    x.SetTitleOffset(1.5)
    x.CenterTitle()
    c.BuildLegend(0.65, 0.7, 0.9, 0.9)
    c.SaveAs('reg1.png')

    hlist = analysisManager.GetProcStack(region='4j2b')
    hs = THStack('j4b1', '>=4 jets, 2 b-tag; H_{T} [GeV]')
    for h in hlist:
        hs.Add(h)
    hs.Draw('hist pfc plc')
    c.Draw()
    x = hs.GetXaxis()
    x.SetTitleOffset(1.5)
    x.CenterTitle()
    c.BuildLegend(0.65, 0.7, 0.9, 0.9)
    c.SaveAs('reg2.png')

    freshstack = analysisManager.GetVarStack(region='4j1b')
    hs = THStack('j4b1btag', 'btag-variations ; H_{T} [GeV]')
    for h in freshstack:
        ptr = h.Rebin(2, h.GetTitle())
        ptr.SetFillColor(0)
        ptr.SetLineWidth(1)
        hs.Add(ptr)
    hs.Draw('hist nostack')
    c.Draw()
    x = hs.GetXaxis()
    x.SetRangeUser(120, 500)
    x.SetTitleOffset(1.5)
    x.CenterTitle()
    c.BuildLegend(0.65, 0.7, 0.9, 0.9)
    c.SaveAs('btag.png')

    output = ROOT.TFile.Open(FILE, 'RECREATE')
    for process in analysisManager:
        for variation in analysisManager[process]:
            for region in analysisManager[process][variation]:
                hist_name = f"{region}_{process}_{variation}" if variation != 'nominal' else f"{region}_{process}"
                hist = analysisManager[process][variation][region]
                if not isinstance(hist, ROOT.TH1D): hist = hist.GetValue() #this this a bag
                if hist.IsZombie(): raise TypeError(hist_name)
                hist_sliced = ROOT.Slice(hist, 120, 550)
                hist_binned = hist_sliced.Rebin(2, hist.GetTitle())
                output.WriteObject(hist_binned, hist_name)
    output.Close()

def main():
    results = analyse()
    make_plots(results)

if __name__ == "__main__":
    raise SystemExit(main())
