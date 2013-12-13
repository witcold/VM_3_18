program ProjectMatrix;

uses
  Forms,
  UnitMatrix in 'UnitMatrix.pas' {FormMain};

{$R *.res}

begin
  Application.Initialize;
  Application.CreateForm(TFormMain, FormMain);
  Application.Run;
end.
