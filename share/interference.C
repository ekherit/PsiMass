{
  gSystem->CompileMacro("interference.cpp","kO");
  TCanvas * c1 = new TCanvas;
  double deg[5] = {10,20,30,40,60};
  TLegend * l1 = new TLegend(0.8,0.6,1.0,1.0);
  for(int i=0;i<5;i++)
  {
    char buf[1024];
    sprintf(buf,"f%d",i);
    TF1 * f = new TF1(buf,&ee_interference_vs_W,3680,3690,1);
    f->SetParameter(0,deg[i]/180.*3.1415926535);
    f->SetLineColor(i+1);
    if(i==0) f->Draw();
    else f->Draw("same");
    sprintf(buf,"#theta=%2.1f",deg[i]);
    l1->AddEntry(f,buf,"l");
  }
  l1->Draw();
  TCanvas * c2 = new TCanvas;
}
