"""In-process Atoms object cache with simple string keys.

Thread-safe LRU map storing ASE Atoms instances for reuse between requests.
"""

from __future__ import annotations

import os
import threading
import uuid
from collections import OrderedDict
from typing import Optional

from fastapi import HTTPException


class _AtomsCache:
    def __init__(self, max_size: int = 512):
        self._lock = threading.Lock()
        self._max_size = max(1, int(max_size))
        self._store: "OrderedDict[str, object]" = OrderedDict()

    def get(self, key: str):
        with self._lock:
            try:
                atoms = self._store.pop(key)
            except KeyError:
                raise HTTPException(status_code=404, detail="cache key not found")
            # move-to-end for LRU recency
            self._store[key] = atoms
            return atoms

    def put(self, atoms, key: Optional[str] = None) -> str:
        """Insert or update an atoms object.

        If key is None, generate a new UUID4 string. Returns the key used.
        """
        if key is None:
            key = uuid.uuid4().hex
        with self._lock:
            if key in self._store:
                # Overwrite and move to end
                self._store.pop(key, None)
            self._store[key] = atoms
            # Enforce capacity
            while len(self._store) > self._max_size:
                self._store.popitem(last=False)
        return key

    def clear(self):
        with self._lock:
            self._store.clear()

    def __len__(self):
        with self._lock:
            return len(self._store)


_DEFAULT_SIZE = int(os.environ.get("UMA_ATOMS_CACHE_SIZE", "512"))
CACHE = _AtomsCache(_DEFAULT_SIZE)


def cache_get(key: str):
    return CACHE.get(key)


def cache_put(atoms, key: Optional[str] = None) -> str:
    return CACHE.put(atoms, key)


def cache_clear():
    CACHE.clear()


__all__ = ["cache_get", "cache_put", "cache_clear", "CACHE"]
