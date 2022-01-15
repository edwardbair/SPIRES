function RasterReprojectionMemoryCheck(neededMemory)

persistent already % send the annoying warning to the Mac/Unix folk just once
if ispc
    [~,sysview] = memory;
    availablePhysicalMemory = sysview.PhysicalMemory.Available;
    availableVirtualMemory = sysview.SystemMemory.Available;
    if neededMemory>availableVirtualMemory
        error('not enough memory available, virtual=%g, needed=%g (so break the problem up)',...
            availableVirtualMemory,neededMemory)
    end
    if neededMemory>availablePhysicalMemory
        warning('not enough physical memory available, physical=%g, needed=%g (using virtual memory may slow processing)',...
            availablePhysicalMemory,neededMemory)
    end
    
elseif ismac
    if isempty(already)
        warning('MATLAB MemoryCheck function not available on a Mac, your code could blow up if input raster or 3D object is large, good luck')
        already = true;
    end
elseif isunix
    if isempty(already)
        warning('MATLAB MemoryCheck function not available on Unix, your code could blow up if input raster or 3D object is large, good luck')
        already = true;
    end
end
end