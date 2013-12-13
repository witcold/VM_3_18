unit UnitMatrix;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, Grids, Buttons;

const
  Ranks: array [1..2] of Integer = (10, 20);
  Ranges: array [1..2] of Integer = (2, 50);
  Epsilons: array [1..3] of Double = (1E-5, 1E-7, 1E-9);
  Tests = 10;

type
  TType = real;
  TMatrix = array of array of TType;
  TVector = array of TType;

  TFormMain = class(TForm)
    StringGridResult: TStringGrid;
    ButtonTest: TButton;
    procedure FormCreate(Sender: TObject);
    procedure ButtonTestClick(Sender: TObject);
    Procedure Jacobi (n:integer; A:TMatrix; max_aij:TType; var k:integer; var sr_oc_toch_L:TType; var sr_mera_toch_r:TType; lyambda:TVector);
    procedure GenerateA(var A:TMatrix; n, min_L, max_L:integer; var y:TVector);

  private
    { Private declarations }
  public
    { Public declarations }
  end;


var
  FormMain: TFormMain;
  A, H: TMatrix;

implementation

uses Math;

{$R *.dfm}

procedure TFormMain.FormCreate(Sender: TObject);
var
  i: Integer;
begin
  with StringGridResult do begin
    Cells[0,0] := '№';
    Cells[1,0] := 'Размерность с-мы N';
    Cells[2,0] := 'Диапазон значений L';
    Cells[3,0] := 'Max Aij';
    Cells[4,0] := 'Ср. число итераций';
    Cells[5,0] := 'Ср. оценка точности L';
    Cells[6,0] := 'Ср. мера точности r';
    for i:=1 to 12 do
      Cells[0,i] := IntToStr(i);
  end;
end;

procedure InitM (var a:TMatrix;n:integer);
var
  i: Integer;
begin
  setlength(a,n);
  for i:=0 to n-1 do
    setlength (a[i], n);
end;

function MultiVect (A,B:TVector;n:integer):TMatrix;
var
  i,j: Integer;
begin
  InitM(Result,n);
  for i:=0 to n-1 do
    for j:=0 to n-1 do
      Result[i,j]:=A[i]*B[j];
end;

function trans (a:TMatrix;n:integer):TMatrix;
var 
  i,j:integer;  
  tr:TMatrix;
begin
  InitM(tr,n);
  for i:=0 to n-1 do
    for j:=0 to n-1 do
      tr[i,j]:=a[j,i];
  Result:=tr;
end;

function MultiMatr (A,B:TMatrix;n:integer):TMatrix;
var
  i,j,k: Integer;
begin
 InitM(Result,n);
 for i:=0 to n-1 do
  for j:=0 to n-1 do
   begin
     Result[i,j]:=0;
     for k:=0 to n-1 do
      Result[i,j] := Result[i,j] + A[i][k]*B[k][j];
   end;
end;

procedure Upor (var M:TVector; n:integer);
var 
  i,j: integer;
  t: TType;
begin
  for j:=0 to n-2 do
    for i:=0 to n-j-2 do
      if M[i] > M[i+1] then begin
        t := M[i];
        M[i] := M[i+1];
        M[i+1] := t
      end;
end;

procedure TFormMain.GenerateA(var A:TMatrix; n, min_L, max_L:integer; var y:TVector);
var
  i,j:integer;
  x:TVector;
  diag: TMatrix;
  sum:TType;
begin
  //генерация собственных значений лямбда
  SetLength(y,n);
  randomize;
  for i:=0 to n-1 do
    y[i]:=random(max_L)-min_L+random;

  Upor (y,n);

  //диагональная матрица лямбда, на ее диагонали стоят собственные значения
  InitM(diag,n);
  for i:=0 to n-1 do
    diag [i,i]:=y[i];

  //генерация и нормирование вектора X
  SetLength(x,n);
  sum:=0;
  for i:=0 to n-1 do begin
    x[i]:=random(5)+1;
    sum:=sum+x[i]*x[i];
  end;
  sum := 1/sqrt(sum);
  for i:=0 to n-1 do
    x[i]:=x[i]*sum;


  //получение матрицы Хаусхолдера
  InitM(H,n);
  for i:=0 to n-1 do
    for j:=0 to n-1 do
     if i=j then
       H[i,j]:=1-2*X[i]*X[j]
     else
       H[i,j]:=-2*X[i]*X[j];
  //получение матрицы А
  A:=MultiMatr(MultiMatr(H,diag,n),trans(H,n),n);
end;

procedure Method (var Bpred:TMatrix; var Bk:TMatrix; var T:TMatrix; n:integer);
var
  imax, jmax, i,j:integer;
  max, p, q,d,r,c,s,pr:TType;
  Tnew:TMatrix;
begin
  imax:=0;
  jmax:=1;
  max:=abs(Bpred[0,1]);
  for i:=0 to n-1 do
   for j:=i+1 to n-1 do
     if abs(Bpred[i,j])>max then begin
       max:=abs(Bpred[i,j]);
       imax:=i;
       jmax:=j
     end;

  p:= 2*BPred[imax,jmax];
  q:=Bpred[imax, imax]-Bpred[jmax,jmax];
  d:=sqrt(sqr(p)+sqr(q));
  if q <> 0 then begin
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

  //проверка
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
  T:=MultiMatr(T, Tnew, n);
end;

procedure Pogresh (lambda, sobst:TVector; var pogr: TType; n: Integer);
var
  i: Integer;
begin
  pogr:=0;
  for i:=0 to n-1 do
    if abs(lambda[i]-sobst[i]) > pogr then
      pogr:=abs(lambda[i] - sobst[i]);
end;

Procedure TFormMain.Jacobi (n:integer; A:TMatrix; max_aij:TType; var k:integer; var sr_oc_toch_L:TType; var sr_mera_toch_r:TType; lyambda:TVector);
//проверка, является ли матрица диагональной
  function DiagMatr (var A:TMatrix; n:integer):boolean;
  var  ok:boolean; i,j:integer;
  begin
    ok:=true;
    i:=0; j:=0;
    while ok and (i<=n-1) do begin
      while ok and (j<=n-1) do begin
        if i=j then ok:=true
        else if abs(A[i,j])<max_aij  then ok:=true
        else ok:=false;
        j:=j+1;
      end;
      i:=i+1;
      j:=0;
    end;
    DiagMatr:=ok;
  end;

var
  i,j: Integer;
  M:integer;
  T:TMatrix;
  sobst:TVector;
  max:TType;
  Bpred, bk, diag, TD,AT:TMatrix;
begin
  M:=10000;
  //первоначально единичная матрица
  InitM (T, n);
  for i:=0 to n-1 do
    T[i,i]:=1;
  SetLength(sobst, n);
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
      if abs(TD[i,j]-AT[i,j])>max then max:=abs(TD[i,j]-AT[i,j]);
  sr_mera_toch_r:=max;

  Upor (sobst, n);
  pogresh (lyambda, sobst,sr_oc_toch_L, n);
end;

procedure TFormMain.ButtonTestClick(Sender: TObject);
var
  i, j, k, t, line, iterations: integer;
  sr_oc_toch_L, sr_mera_toch_r: TType;
  lambda: TVector;
  matr_oc_L, matr_mera: TType;
  matr_k: integer;
begin
  line := 1;
  for i := 1 to Length(Ranks) do
    for j := 1 to Length(Ranges) do
      for k := 1 to Length(Epsilons) do begin
        iterations := 0;
        sr_oc_toch_L:=0;
        sr_mera_toch_r:=0;
        for t := 1 to Tests do begin
          matr_k:=0;
          matr_oc_L:=0;
          matr_mera:=0;
          InitM(A, Ranks[i]);
          GenerateA(A, Ranks[i], -Ranges[j], Ranges[j], lambda);
          Jacobi(Ranks[i], A, Epsilons[k],matr_k,matr_oc_L, matr_mera, lambda);
          sr_oc_toch_L:=sr_oc_toch_L+matr_oc_L;
          sr_mera_toch_r:=sr_mera_toch_r+matr_mera;
          iterations := iterations + matr_k;
        end;
        sr_oc_toch_L:=sr_oc_toch_L/Tests;
        sr_mera_toch_r:=sr_mera_toch_r/Tests;
        iterations := round(iterations / Tests);
        with StringGridResult do begin
          Cells[1, line] := IntToStr(Ranks[i]);
          Cells[2, line] := IntToStr(-Ranges[j])+'...'+IntToStr(Ranges[j]);
          Cells[3, line] := FloatToStr(Epsilons[k]);
          Cells[4, line] := IntToStr(iterations);
          Cells[5, line] := FloatToStr(sr_oc_toch_L);
          Cells[6, line] := FloatToStr(sr_mera_toch_r);
        end;
        Application.ProcessMessages;
        line := line + 1;
      end;
end;

end.
