function LoadTSeriestoBruker(TSeriesENVFile)

pl = actxserver('PrairieLink.Application');
pl.Connect();
pl.SendScriptCommands(['-TSeriesLoad ' TSeriesENVFile.folder '\' TSeriesENVFile.name]);
pl.Disconnect();
delete(pl);
clear pl