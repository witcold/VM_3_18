unit UnitMatrix;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, Grids, Buttons;

const
  Ranks: array [1..2] of Integer = (10, 20);
  Ranges: array [1..2] of Integer = (2, 50);
  Epsilons: array [1..3] of Double = (1E-5, 1E-7, 1E-9);

type
  TElem = real;
  TMatr = array of array of TElem;
  TVector = array of TElem;
  TFormMain = class(TForm)
    StringGridResult: TStringGrid;
    ButtonTest: TButton;
  procedure FormCreate(Sender: TObject);
  procedure ButtonTestClick(Sender: TObject);
  Procedure Jacobi (n:integer; A:TMAtr; max_aij:TElem; var k:integer; var sr_oc_toch_L:TElem; var sr_mera_toch_r:TElem; lyambda:TVector);
  procedure GenerateA(var A:Tmatr; n, min_L, max_L:integer; var y:TVector);

  private
    { Private declarations }
  public
    { Public declarations }
  end;


var
  FormMain: TFormMain;
  A, H: TMatr;

implementation

uses Math;

{$R *.dfm}

procedure TFormMain.FormCreate(Sender: TObject);
var
  i,j: integer;
begin
  StringGridResult.Cells[0,0] := '¹';
  StringGridResult.Cells[1,0] := 'Ðàçìåðíîñòü ñ-ìû N';
  StringGridResult.Cells[2,0] := 'Äèàïàçîí çíà÷åíèé L';
  StringGridResult.Cells[3,0] := '     Max Aij';
  StringGridResult.Cells[4,0] := 'Ñð. ÷èñëî èòåðàöèé';
  StringGridResult.Cells[5,0] := 'Ñð. îöåíêà òî÷íîñòè L';
  StringGridResult.Cells[6,0] := 'Ñð. ìåðà òî÷íîñòè r';
  for i:=1 to 12 do
    StringGridResult.Cells[0,i] := IntToStr(i);
end;

procedure InitV (var v:TVector; n:integer);
var
  i:integer;
begin
  setlength (v, n);
  for i:=0 to n-1 do
    v[i]:=0;
end;

procedure InitM (var a:TMatr;n:integer);
var
  i,j:integer;
begin
  setlength(a,n);
  for i:=0 to n-1 do begin
    setlength (a[i], n);
  for j:=0 to n-1 do
    a[i][j]:=0;
  end;
end;

function MultiVect (A,B:TVector;n:integer):TMatr;
var
  i,j,k:integer;
  sum:real;
  c:TMatr;
begin
  InitM(C,n);
  for i:=0 to n-1 do
    for j:=0 to n-1 do
      C[i,j]:=A[i]*B[j];
  MultiVect:=C;
end;

function trans (a:TMatr;n:integer):TMatr;
var
  i,j:integer;
  tr:TMatr;
begin
   InitM(tr,n);
   for i:=0 to n-1 do
    for j:=0 to n-1 do
     tr[i,j]:=a[j,i];
   trans:=tr;
end;

function MultiMatr (A,B:TMatr;n:integer):TMatr;
var
  i,j,k:integer;
  sum:real;
  c:TMatr;
begin
  InitM(C,n);
  for i:=0 to n-1 do
    for j:=0 to n-1 do begin
      sum:=0;
      for k:=0 to n-1 do
        sum:=sum+A[i,k]*B[k,j];
      C[i][j]:=sum;
    end;
  MultiMatr:=C;
end;

//ñîðòèðîâêà âåêòîðà (çà÷åì?)
procedure Upor (var M:TVector; n:integer);
var
  i,j:integer;
  t:TElem;
begin
  for j := 0 to n-2 do
   for i := 0 to n-j-2 do
     if M[i] > M[i+1] then begin
       t := M[i];
       M[i] := M[i+1];
       M[i+1] := t
     end;
end;

procedure TFormMain.GenerateA(var A:Tmatr; n, min_L, max_L:integer; var y:TVector);
var
    i,j:integer;
    x:TVector;
    diag,XXt, E, T:TMatr;
    sum:TElem;
begin
  //ãåíåðàöèÿ ñîáñòâåííûõ çíà÷åíèé ëÿìáäà
  InitV(y,n);
  randomize;
  for i:=0 to n-1 do
    y[i]:=random(max_L)-min_L+random;
  Upor (y,n);

  //äèàãîíàëüíàÿ ìàòðèöà ëÿìáäà, íà åå äèàãîíàëè ñòîÿò ñîáñòâåííûå çíà÷åíèÿ
  InitM(diag,n);
  for i:=0 to n-1 do
    diag [i,i]:=y[i];

  //ãåíåðàöèÿ âåêòîðà
  InitV(x,n);
  for i:=0 to n-1 do
    x[i]:=random(5)+1;

  //íîðìèðîâàíèå âåêòîðà x
  sum := 0;
  for i:=0 to n-1 do
    sum := sum+x[i]*x[i];
  for i:=0 to n-1 do
    x[i]:=x[i]/sqrt(sum);

  //ïåðåìíîæåíèå âåêòîðà õ íà òðàíñïîíèðîâàííûé õ
  InitM(XXt,n);
  XXt:=MultiVect(X, X,n);

  //åäèíè÷íàÿ ìàòðèöà
  InitM(E,n);
  for i:=0 to n-1 do
   E[i,i]:=1;

  //ïîëó÷åíèå ìàòðèöû Õàóñõîëäåðà
  InitM(H,n);
  for i:=0 to n-1 do
    for j:=0 to n-1 do
      H[i,j]:= E[i,j]-2*XXt[i,j]; //èçìåíèòü íà ñîîòâåòñòâóþùèé ýëåìåíò

  //ïîëó÷åíèå ìàòðèöû À
  A:=MultiMatr(MultiMatr(H,diag,n),trans(H,n),n);
end;



procedure Method (var Bpred:TMatr; var Bk:TMatr; var T:TMatr; n:integer);
var
  imax, jmax, i,j:integer;
  max, p, q,d,r,c,s,pr:TElem;
  cnew,snew:TElem;
  Tnew:TMatr;
begin
  imax:=0;
  jmax:=1;
  max:=abs(Bpred[0,1]);
  for i:=0 to n-1 do
    for j:=i+1 to n-1 do
      if abs(Bpred[i,j])>max then begin
        max := abs(Bpred[i,j]);
        imax := i;
        jmax := j
      end;
  p:= 2*BPred[imax,jmax];
  q:=Bpred[imax, imax]-Bpred[jmax,jmax];
  d:=sqrt(sqr(p)+sqr(q));
  if q<>0 then begin
    r:=abs(q)/(2*d);
    c:=sqrt(0.5+r);
    s:=sqrt(0.5-r)*sign(p*q)
  end
  else begin
    c:=sqrt(2)/2;
    s:=c;
  end;
  Bk[imax,imax]:=sqr(c)*Bpred[imax,imax]+sqr(s)*Bpred[jmax,jmax]+2*c*s*Bpred[imax,jmax];
  Bk[jmax,jmax]:=sqr(s)*Bpred[imax,imax]+sqr(c)*Bpred[jmax,jmax]-2*c*s*Bpred[imax,jmax];
  Bk[imax,jmax]:=0;
  Bk[jmax,imax]:=0;
  //ïðîâåðêà
  pr:=(sqr(c)-sqr(s))*Bpred[imax, jmax]+c*s*(Bpred[jmax, imax]-Bpred[imax,imax]);
  for i:=0 to n-1 do
    if (i<>imax) and (i<>jmax) then begin
      Bk[imax, i]:=c*Bpred[i,imax]+s*Bpred[i,jmax];
      Bk[i,imax]:=Bk[imax,i];
      Bk[jmax, i]:=-s*Bpred[i,imax]+c*Bpred[i,jmax];
      Bk[i,jmax]:= Bk[jmax, i];
    end;
  for i:=0 to n-1 do
    for j:=0 to n-1 do
      if (i<>imax) and (j<>jmax)and (i<>jmax)and (j<>imax) then
        Bk[i,j]:=Bpred[i,j];
  InitM (Tnew, n);
  for i:=0 to n-1 do
    Tnew[i,i]:=1;
  Tnew[imax, imax]:=c;
  Tnew[imax, jmax]:=-s;
  Tnew[jmax, jmax]:=c;
  Tnew[jmax, imax]:=s;
  T:=MultiMatr (T, Tnew, n);
end;

procedure pogresh (lyambda, sobst:TVector; var pogr:TElem; n:integer);
var
  i,j,js:integer;
  find:array of boolean;
  min:TVector;
begin
  pogr := 0;
  for i := 0 to n-1 do
   if abs(lyambda[i]-sobst[i])>pogr then
     pogr:=abs(lyambda[i]-sobst[i]);
end;

Procedure TFormMain.Jacobi (n:integer; A:TMAtr; max_aij:TElem; var k:integer; var sr_oc_toch_L:TElem; var sr_mera_toch_r:TElem; lyambda:TVector);
//ïðîâåðêà, ÿâëÿåòñÿ ëè ìàòðèöà äèàãîíàëüíîé
  function DiagMatr (var A:Tmatr; n:integer):boolean;
  var
    ok:boolean;
    i,j:integer;
  begin
    ok:=true;
    i:=0;
    j:=0;
    while ok and (i<=n-1) do begin
      while ok and (j<=n-1) do begin
        if i=j then
          ok:=true
        else
          if abs(A[i,j])<max_aij  then
            ok:=true
          else ok:=false;
        j:=j+1;
      end;
      i:=i+1;
      j:=0;
    end;
    DiagMatr:=ok;
  end;

var
  i,j,w:integer;
  M:integer;
  T:Tmatr;
  sobst:TVector;
  r:TElem;
  max:TElem;
  Bpred, bk:TMatr;
  diag, TD,AT:TMatr;
begin
  M:=10000;
  //ïåðâîíà÷àëüíî åäèíè÷íàÿ ìàòðèöà
  InitM (T, n);
  for i:=0 to n-1 do
    T[i,i]:=1;
  InitV(sobst, n);
  r:=0;
  InitM(Bpred, n);
  InitM(Bk, n);
  for i:=0 to n-1 do
    for j:=0 to n-1 do begin
      Bpred[i,j]:=A[i,j];
      Bk[i,j]:=A[i,j];
    end;
  while not DiagMatr (Bpred, n) and (k<=M) do begin
    Method (BPred, Bk, T, n);
    k:=k+1;
    for i:=0 to n-1 do
      for j:=0 to n-1 do
        Bpred[i,j]:=Bk[i,j];
  end;
  InitM(diag, n);
  for i:=0 to n-1 do begin
    sobst[i]:=bk[i,i];
    diag[i,i]:=bk[i,i];
  end;
  InitM(TD, n);
  InitM(AT, n);
  TD:=MultiMatr(T,diag, n);
  AT:=MultiMatr(A,T, n);
  max:=abs(TD[0,0]-AT[0,0]);
  for i:=0 to n-1 do
    for j:=0 to n-1 do
      if abs(TD[i,j]-AT[i,j])>max then
        max:=abs(TD[i,j]-AT[i,j]);
  sr_mera_toch_r:=max;
  Upor (sobst, n);
  pogresh (lyambda, sobst,sr_oc_toch_L, n);
end;

procedure TFormMain.ButtonTestClick(Sender: TObject);
var
  k, i: integer;
  pogr, sr_oc_toch_L, sr_mera_toch_r: telem;
  lyambda: TVector;
  matr_oc_L, matr_mera:array[1..10] of TElem;
  matr_k:array[1..10] of integer;
  rank, range, eps, line: Integer;
begin
  line := 1;
  for rank := 1 to Length(Ranks) do
    for range := 1 to Length(Ranges) do
      for eps := 1 to Length(Epsilons) do begin
        k:=0;
        sr_oc_toch_L:=0;
        sr_mera_toch_r:=0;
        for i:=1 to 10 do begin
          matr_k[i]:=0;
          matr_oc_L[i]:=0;
          matr_mera[i]:=0;
        end;
        for i:=1 to 10 do begin
          InitM(A, Ranks[rank]);
          GenerateA(A, Ranks[rank], -Ranges[range], Ranges[range], lyambda);
          Jacobi(Ranks[rank], A, Epsilons[eps],matr_k[i],matr_oc_L[i], matr_mera[i], lyambda);
        end;
        for i:=1 to 10 do begin
          sr_oc_toch_L:=sr_oc_toch_L+matr_oc_L[i];
          sr_mera_toch_r:=sr_mera_toch_r+matr_mera[i];
          k:=k+matr_k[i];
        end;
        sr_oc_toch_L:=sr_oc_toch_L/10;
        sr_mera_toch_r:=sr_mera_toch_r/10;
        k:=round (k/10);
        StringGridResult.Cells[1, line] := IntToStr(Ranks[rank]);
        StringGridResult.Cells[2, line] := IntToStr(-Ranges[range])+'...'+IntToStr(Ranges[range]);
        StringGridResult.Cells[3, line] := FloatToStr(Epsilons[eps]);
        StringGridResult.Cells[4, line] := IntToStr(k);
        StringGridResult.Cells[5, line] := FloatToStr(sr_oc_toch_L);
        StringGridResult.Cells[6, line] := FloatToStr(sr_mera_toch_r);
        line := line + 1;
      end;
end;

end.
