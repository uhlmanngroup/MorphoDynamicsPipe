from __future__ import annotations
import torch


def get_size_chunk(device_mem: int, chunk_size: int):
    """get the number of `chunk_size` possible to fit inside the `device_mem`

    Parameters
    ----------
    device_mem : int
        device memory size
    chunk_size : int, optional
        size of the chunk we want, in GiB
        This parameter is dependant on the model to use

    Notes
    -----
    .. warning:: Size
        Size of the chunks are in GiB no GB (decimal bytes)


    Returns
    -------
    int
        number of chunks presents in this device_mem. `int()` floors it to lower integer
    """
    chunk_size = chunk_size * (1000**3)
    return int(device_mem / chunk_size)


def get_gpu_chunks(chunk_size: int = 4):
    """Parse the gpu devices and count the chunks of memory possible to run
    on the GPU VRAM

    Parameters
    ----------
    chunk_size : int, optional
        size of the chunk we want, in GiB
        This parameter is dependant on the model to use

    Returns
    -------
    int
        number of chunks of memory available across the gpu(s) vram
    Notes
    -----
    .. warning:: Size
        Size of the chunks are in GiB no GB (decimal bytes)
    """
    if not torch.cuda.is_available() or torch.cuda.device_count() == 0:
        return 0
    count_ = 0
    for device_index in range(torch.cuda.device_count()):
        if device_index is None or device_index == -1:
            print(f"Error on cuda device: {device_index}")
            continue
        current_device = torch.cuda.device(device_index)
        count_ += get_size_chunk(
            torch.cuda.get_device_properties(current_device).total_memory,
            chunk_size=chunk_size,
        )
    return count_


if __name__ == "__main__":
    print(f"available cellpose threads = {get_gpu_chunks()}")
